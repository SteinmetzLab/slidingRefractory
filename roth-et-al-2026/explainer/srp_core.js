/*
 * srp_core.js -- numerical core of the Sliding RP explorer.
 *
 * A JavaScript re-expression of the parts of the slidingRP package the explorer
 * needs, for a single unit:
 *   - the manuscript's spike-train model (hard-RP renewal base neuron plus a
 *     Poisson contaminating source, contamination = fraction of the total);
 *   - computeACG (histdiff-equivalent, 1/30000 s bins, 0-10 ms);
 *   - the Sliding RP confidence row, acceptance, first/last accepted tau_r and
 *     C_min, exactly as slidingRP.metrics.slidingRP computes them;
 *   - the confidence matrix, as computeMatrix;
 *   - Hill-Llobet at a fixed RP, as matlab/RPmetric_Classic.m ('Llobet').
 *
 * Checked against the Python package by test_core.js (run: node test_core.js
 * after python export_reference.py). Works in a browser (window.SRP) and in node
 * (module.exports).
 */
(function (root) {
  'use strict';

  const SAMPLE_RATE = 30000;
  const BIN = 1 / SAMPLE_RATE;
  const N_BINS = 300;                        // 0-10 ms
  const MAX_LAG = BIN * N_BINS;
  const RP_CENTERS = new Float64Array(N_BINS);
  const REF_DUR = new Float64Array(N_BINS);
  for (let k = 0; k < N_BINS; k++) {
    RP_CENTERS[k] = k * BIN + BIN / 2;       // numpy: arange(0, 10ms, bin) + bin/2
    REF_DUR[k] = RP_CENTERS[k] + BIN / 2;    // window duration used for V_e
  }
  const CONT_GRID = [];
  for (let i = 1; i <= 70; i++) CONT_GRID.push(0.5 * i);   // 0.5:0.5:35 (%)

  // ---------------------------------------------------------------- RNG ----
  // sfc32 seeded through splitmix32; doubles with 53 random bits.
  function makeRng(seed) {
    let s = seed >>> 0;
    function splitmix() {
      s = (s + 0x9e3779b9) >>> 0;
      let z = s;
      z = Math.imul(z ^ (z >>> 16), 0x85ebca6b) >>> 0;
      z = Math.imul(z ^ (z >>> 13), 0xc2b2ae35) >>> 0;
      return (z ^ (z >>> 16)) >>> 0;
    }
    let a = splitmix(), b = splitmix(), c = splitmix(), d = splitmix();
    function u32() {
      a >>>= 0; b >>>= 0; c >>>= 0; d >>>= 0;
      let t = (a + b) | 0;
      a = b ^ (b >>> 9);
      b = (c + (c << 3)) | 0;
      c = (c << 21) | (c >>> 11);
      d = (d + 1) | 0;
      t = (t + d) | 0;
      c = (c + t) | 0;
      return t >>> 0;
    }
    for (let i = 0; i < 12; i++) u32();
    return function () {                     // uniform on [0, 1)
      return (u32() * 2097152 + (u32() >>> 11)) / 9007199254740992;
    };
  }

  // ------------------------------------------------------- spike trains ----
  /** Renewal process ISI = rp + Exp(mu), rate-corrected for the dead time
   *  (the manuscript's genST). Returns sorted times in [0, duration). */
  function genHardRP(rate, duration, rp, rng) {
    if (!(rate > 0)) return new Float64Array(0);
    if (rp * rate >= 1) throw new Error('rate too high for this refractory period');
    const mu = (1 - rp * rate) / rate;
    let cap = Math.max(16, Math.ceil(rate * duration * 1.2 + 10 * Math.sqrt(rate * duration + 1)));
    let out = new Float64Array(cap);
    let n = 0, t = 0;
    for (;;) {
      t += rp - mu * Math.log(1 - rng());
      if (t >= duration) break;
      if (n === cap) {
        const bigger = new Float64Array(cap * 2);
        bigger.set(out);
        out = bigger;
        cap *= 2;
      }
      out[n++] = t;
    }
    return out.subarray(0, n);
  }

  function mergeSorted(a, b) {
    const out = new Float64Array(a.length + b.length);
    let i = 0, j = 0, k = 0;
    while (i < a.length && j < b.length) out[k++] = a[i] <= b[j] ? a[i++] : b[j++];
    while (i < a.length) out[k++] = a[i++];
    while (j < b.length) out[k++] = b[j++];
    return out;
  }

  /** The manuscript's standard model: contProp is the fraction of the total
   *  train contributed by a Poisson source with no refractory period. */
  function makeTrain(totalRate, contProp, duration, rp, rng) {
    const base = genHardRP((1 - contProp) * totalRate, duration, rp, rng);
    const cont = genHardRP(contProp * totalRate, duration, 0, rng);
    return { st: mergeSorted(base, cont), nBase: base.length, nCont: cont.length };
  }

  // ---------------------------------------------------------------- ACG ----
  /** histdiff-equivalent: bin = floor(dt / BIN) for ordered pairs with
   *  0 < dt < 10 ms. Spike times must be sorted. */
  function computeACG(st) {
    const acg = new Float64Array(N_BINS);
    const n = st.length;
    for (let i = 0; i < n; i++) {
      const ti = st[i];
      for (let j = i + 1; j < n; j++) {
        const dt = st[j] - ti;
        if (dt >= MAX_LAG) break;
        if (dt > 0) {
          const k = Math.floor(dt / BIN);
          if (k < N_BINS) acg[k] += 1;
        }
      }
    }
    return acg;
  }

  // ---------------------------------------------------------- Poisson ------
  const LANCZOS = [57.1562356658629235, -59.5979603554754912, 14.1360979747417471,
    -0.491913816097620199, 0.339946499848118887e-4, 0.465236289270485756e-4,
    -0.983744753048795646e-4, 0.158088703224912494e-3, -0.210264441724104883e-3,
    0.217439618115212643e-3, -0.164318106536763890e-3, 0.844182239838527433e-4,
    -0.261908384015814087e-4, 0.368991826595316234e-5];

  function lnGamma(xx) {                     // Numerical Recipes 3rd ed. gammln
    let x = xx, y = xx;
    let tmp = x + 5.24218750000000000;
    tmp = (x + 0.5) * Math.log(tmp) - tmp;
    let ser = 0.999999999999997092;
    for (let j = 0; j < 14; j++) ser += LANCZOS[j] / ++y;
    return tmp + Math.log(2.5066282746310005 * ser / x);
  }

  const ITMAX = 1000000, EPS = 1e-15, FPMIN = 1e-300;

  function gser(a, x) {                      // lower regularized P(a, x), series
    let ap = a, sum = 1 / a, del = sum;
    for (let n = 0; n < ITMAX; n++) {
      ap += 1;
      del *= x / ap;
      sum += del;
      if (Math.abs(del) < Math.abs(sum) * EPS) break;
    }
    return sum * Math.exp(-x + a * Math.log(x) - lnGamma(a));
  }

  function gcf(a, x) {                       // upper regularized Q(a, x), Lentz
    let b = x + 1 - a, c = 1 / FPMIN, d = 1 / b, h = d;
    for (let i = 1; i <= ITMAX; i++) {
      const an = -i * (i - a);
      b += 2;
      d = an * d + b;
      if (Math.abs(d) < FPMIN) d = FPMIN;
      c = b + an / c;
      if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < EPS) break;
    }
    return Math.exp(-x + a * Math.log(x) - lnGamma(a)) * h;
  }

  function gammp(a, x) {
    if (x <= 0) return 0;
    return x < a + 1 ? gser(a, x) : 1 - gcf(a, x);
  }

  /** P(X > k) for X ~ Poisson(lam): the confidence (as a fraction) that the
   *  true contamination is below the level that produced lam. */
  function poissonSF(k, lam) {
    if (!(lam > 0)) return 0;
    return gammp(k + 1, lam);
  }

  function poissonPMF(k, lam) {
    if (!(lam > 0)) return k === 0 ? 1 : 0;
    return Math.exp(k * Math.log(lam) - lam - lnGamma(k + 1));
  }

  /** P(X <= k) for X ~ Poisson(lam). */
  function poissonCDF(k, lam) {
    if (!(lam > 0)) return 1;
    const a = k + 1;
    return lam < a + 1 ? 1 - gser(a, lam) : gcf(a, lam);
  }

  // ------------------------------------------------------ the metric -------
  function cumsum(acg) {
    const out = new Float64Array(acg.length);
    let s = 0;
    for (let k = 0; k < acg.length; k++) { s += acg[k]; out[k] = s; }
    return out;
  }

  /** Llobet expected violations for every window, at contamination cPct. */
  function expectedViolations(n, recDur, cPct) {
    const c = cPct / 100;
    const nc = c * n, nb = (1 - c) * n;
    const out = new Float64Array(N_BINS);
    for (let k = 0; k < N_BINS; k++) out[k] = 2 * REF_DUR[k] / recDur * nc * (nb + (nc - 1) / 2);
    return out;
  }

  /** Confidence (%) at every window for one contamination level. */
  function confidenceRow(obs, n, recDur, cPct) {
    const lam = expectedViolations(n, recDur, cPct);
    const out = new Float64Array(N_BINS);
    for (let k = 0; k < N_BINS; k++) out[k] = 100 * poissonSF(obs[k], lam[k]);
    return out;
  }

  function confidenceMatrix(obs, n, recDur, grid) {
    return (grid || CONT_GRID).map(function (c) { return confidenceRow(obs, n, recDur, c); });
  }

  /** Solve P(X > k; lam) = g for lam (the chi-square quantile route in the
   *  package, done by bisection so it needs only the CDF above). */
  function lambdaForConfidence(k, g) {
    let lo = 0, hi = k + 1 + 10 * Math.sqrt(k + 1) + 10;
    while (poissonSF(k, hi) < g) hi *= 2;
    for (let it = 0; it < 200; it++) {
      const mid = 0.5 * (lo + hi);
      if (poissonSF(k, mid) < g) lo = mid; else hi = mid;
      if (hi - lo <= 1e-14 * hi) break;
    }
    return 0.5 * (lo + hi);
  }

  /** Everything slidingRP() reports for one unit, from its ACG. */
  function slidingRP(acg, n, recDur, opts) {
    const o = Object.assign({ contThresh: 10, confThresh: 90, tauMin: 0.0005,
      withCmin: true }, opts || {});
    const obs = cumsum(acg);
    const row = confidenceRow(obs, n, recDur, o.contThresh);
    let maxConf = -Infinity, kMax = -1, kFirst = -1, kLast = -1;
    for (let k = 0; k < N_BINS; k++) {
      if (!(RP_CENTERS[k] > o.tauMin)) continue;
      if (row[k] > maxConf) { maxConf = row[k]; kMax = k; }
      if (row[k] >= o.confThresh) { if (kFirst < 0) kFirst = k; kLast = k; }
    }
    if (kMax < 0) maxConf = 0;
    const res = { obs: obs, row: row, maxConf: maxConf, kMax: kMax,
      passes: maxConf >= o.confThresh, kFirst: kFirst, kLast: kLast,
      cMin: NaN, kCmin: -1 };
    if (o.withCmin) {
      const g = o.confThresh / 100;
      let best = Infinity, kb = -1;
      for (let k = 0; k < N_BINS; k++) {
        if (!(RP_CENTERS[k] > o.tauMin)) continue;
        const lam = lambdaForConfidence(obs[k], g);
        const disc = (n - 0.5) * (n - 0.5) - lam * recDur / REF_DUR[k];
        if (disc < 0) continue;
        const c = ((n - 0.5) - Math.sqrt(disc)) / n * 100;
        if (c < best) { best = c; kb = k; }
      }
      if (kb >= 0 && best <= 35) { res.cMin = best; res.kCmin = kb; }
    }
    return res;
  }

  /** Hill-Llobet point estimate at a fixed RP (seconds). Inclusive bin
   *  selection sum(nACG(1:find(rp > RPdur, 1))) as in RPmetric_Classic.m. */
  function hillLlobet(acg, n, recDur, rpDur, contThresh) {
    let idx = 0;
    while (idx < N_BINS && !(RP_CENTERS[idx] > rpDur)) idx++;
    if (idx >= N_BINS) idx = 0;              // numpy argmax of all-False is 0
    let obs = 0;
    for (let k = 0; k <= idx; k++) obs += acg[k];
    const c = contThresh / 100;
    const nc = n * c, nb = n * (1 - c);
    const expected = 2 * rpDur / recDur * nc * (nb + (nc - 1) / 2);
    const inner = 1 - obs * recDur / (n * n * rpDur);
    const est = inner >= 0 ? 1 - Math.sqrt(inner) : NaN;
    return { obs: obs, expected: expected, est: est, passes: obs <= expected,
      window: RP_CENTERS[idx] + BIN / 2, idx: idx };
  }

  const api = { SAMPLE_RATE: SAMPLE_RATE, BIN: BIN, N_BINS: N_BINS,
    RP_CENTERS: RP_CENTERS, REF_DUR: REF_DUR, CONT_GRID: CONT_GRID,
    makeRng: makeRng, genHardRP: genHardRP, makeTrain: makeTrain,
    mergeSorted: mergeSorted, computeACG: computeACG, lnGamma: lnGamma,
    gammp: gammp, poissonSF: poissonSF, poissonCDF: poissonCDF,
    poissonPMF: poissonPMF, cumsum: cumsum, expectedViolations: expectedViolations,
    confidenceRow: confidenceRow, confidenceMatrix: confidenceMatrix,
    lambdaForConfidence: lambdaForConfidence, slidingRP: slidingRP,
    hillLlobet: hillLlobet };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.SRP = api;
})(typeof window !== 'undefined' ? window : this);
