/*
 * test_core.js -- check srp_core.js against the Python package.
 *
 * Run:  python export_reference.py && node test_core.js
 * Exits non-zero on any failure.
 */
'use strict';
const fs = require('fs');
const path = require('path');
const SRP = require('./srp_core.js');

const REF = 'D:/temp/slidingRP_resub/explainer_ref';
const ref = JSON.parse(fs.readFileSync(path.join(REF, 'reference.json'), 'utf8'));
let failures = 0;

function check(cond, msg) {
  if (!cond) { failures++; console.log('  FAIL  ' + msg); }
}
function close(a, b, tol) {
  if (a === null || b === null || Number.isNaN(a) || Number.isNaN(b)) {
    return (a === null || Number.isNaN(a)) && (b === null || Number.isNaN(b));
  }
  return Math.abs(a - b) <= tol;
}
function readF64(file) {
  const buf = fs.readFileSync(file);
  return new Float64Array(buf.buffer, buf.byteOffset, buf.length / 8);
}

// ---- Poisson tails -----------------------------------------------------
let worstSF = 0, worstCDF = 0;
for (const g of ref.poisson) {
  const sf = SRP.poissonSF(g.k, g.lam), cdf = SRP.poissonCDF(g.k, g.lam);
  worstSF = Math.max(worstSF, Math.abs(sf - g.sf));
  worstCDF = Math.max(worstCDF, Math.abs(cdf - g.cdf));
}
console.log(`Poisson tails: ${ref.poisson.length} points, worst |sf| error ` +
            `${worstSF.toExponential(2)}, worst |cdf| error ${worstCDF.toExponential(2)}`);
check(worstSF < 1e-9 && worstCDF < 1e-9, 'Poisson tail accuracy');

let worstInv = 0;
for (const v of ref.inverse) {
  const lam = SRP.lambdaForConfidence(v.k, v.g);
  worstInv = Math.max(worstInv, Math.abs(lam - v.lam) / v.lam);
}
console.log(`Lambda-for-confidence: worst relative error ${worstInv.toExponential(2)}`);
check(worstInv < 1e-9, 'lambda inverse');

// ---- per-unit cases ----------------------------------------------------
for (const c of ref.cases) {
  const st = readF64(path.join(REF, c.name + '.f64'));
  const acg = SRP.computeACG(st);
  let acgDiff = 0;
  for (let k = 0; k < SRP.N_BINS; k++) acgDiff += Math.abs(acg[k] - c.acg[k]);
  check(acgDiff === 0, `${c.name}: ACG differs in ${acgDiff} counts`);

  const d = SRP.slidingRP(acg, c.n, c.dur);
  check(close(d.maxConf, c.default.max_conf, 1e-9),
        `${c.name}: max conf ${d.maxConf} vs ${c.default.max_conf}`);
  check(d.passes === c.default.passes, `${c.name}: passes`);
  check(close(d.cMin, c.default.min_cont, 1e-7),
        `${c.name}: C_min ${d.cMin} vs ${c.default.min_cont}`);
  const rpMin = d.kCmin >= 0 ? SRP.RP_CENTERS[d.kCmin] : null;
  check(close(rpMin, c.default.rp_min_val, 1e-12),
        `${c.name}: rp at C_min ${rpMin} vs ${c.default.rp_min_val}`);
  const first = d.kFirst >= 0 ? SRP.RP_CENTERS[d.kFirst] : null;
  const last = d.kLast >= 0 ? SRP.RP_CENTERS[d.kLast] : null;
  check(close(first, c.tau_first, 1e-12), `${c.name}: first accepted tau`);
  check(close(last, c.tau_last, 1e-12), `${c.name}: last accepted tau`);

  const a = SRP.slidingRP(acg, c.n, c.dur, { contThresh: 15, confThresh: 80, tauMin: 0.001 });
  check(close(a.maxConf, c.alt.max_conf, 1e-9), `${c.name}: alt max conf`);
  check(a.passes === c.alt.passes, `${c.name}: alt passes`);
  check(close(a.cMin, c.alt.min_cont, 1e-7), `${c.name}: alt C_min ${a.cMin} vs ${c.alt.min_cont}`);
  const aFirst = a.kFirst >= 0 ? SRP.RP_CENTERS[a.kFirst] : null;
  check(close(aFirst, c.alt.tau_first_pass, 1e-12), `${c.name}: alt first accepted tau`);

  const M = SRP.confidenceMatrix(d.obs, c.n, c.dur);
  let worst = 0;
  for (let i = 0; i < M.length; i++)
    for (let k = 0; k < SRP.N_BINS; k++) worst = Math.max(worst, Math.abs(M[i][k] - c.matrix[i][k]));
  check(worst < 1e-7, `${c.name}: matrix max |diff| ${worst}`);

  for (const key of Object.keys(c.hl)) {
    const [rpMs, cont] = key.split('_').map(Number);
    const h = SRP.hillLlobet(acg, c.n, c.dur, rpMs / 1000, cont);
    const r = c.hl[key];
    check(h.passes === r.passes, `${c.name}: HL ${key} passes`);
    check(h.obs === r.obs, `${c.name}: HL ${key} obs ${h.obs} vs ${r.obs}`);
    check(close(h.est, r.est, 1e-12), `${c.name}: HL ${key} estimate`);
  }
  console.log(`${c.name.padEnd(10)} n=${String(c.n).padStart(7)}  max conf ` +
              `${d.maxConf.toFixed(4)}  C_min ${Number.isNaN(d.cMin) ? 'NaN' : d.cMin.toFixed(4)}` +
              `  matrix |diff| ${worst.toExponential(1)}`);
}

// ---- generator, end to end ---------------------------------------------
const s = ref.sim;
const rng = SRP.makeRng(12345);
let nS = 0, nH = 0, rateSum = 0, contSum = 0, minIsiOk = true;
for (let i = 0; i < s.n_trains; i++) {
  const base = SRP.genHardRP((1 - s.cont) * s.rate, s.dur, s.rp, rng);
  for (let j = 1; j < base.length; j++) if (base[j] - base[j - 1] < s.rp - 1e-12) minIsiOk = false;
  const cont = SRP.genHardRP(s.cont * s.rate, s.dur, 0, rng);
  const st = SRP.mergeSorted(base, cont);
  rateSum += st.length / s.dur;
  contSum += cont.length / st.length;
  const acg = SRP.computeACG(st);
  nS += SRP.slidingRP(acg, st.length, s.dur, { withCmin: false }).passes ? 1 : 0;
  nH += SRP.hillLlobet(acg, st.length, s.dur, 0.002, 10).passes ? 1 : 0;
}
const pS = nS / s.n_trains, pH = nH / s.n_trains;
function z(p1, p2, n) {
  const p = (p1 + p2) / 2;
  return (p1 - p2) / Math.sqrt(2 * p * (1 - p) / n + 1e-12);
}
console.log(`Generator: mean rate ${(rateSum / s.n_trains).toFixed(3)} (want ${s.rate}), ` +
            `mean contamination ${(contSum / s.n_trains).toFixed(4)} (want ${s.cont})`);
console.log(`Acceptance over ${s.n_trains} trains: Sliding RP JS ${pS.toFixed(3)} vs ` +
            `Python ${s.p_sliding.toFixed(3)} (z = ${z(pS, s.p_sliding, s.n_trains).toFixed(2)}); ` +
            `HL 2 ms JS ${pH.toFixed(3)} vs Python ${s.p_hl2.toFixed(3)} ` +
            `(z = ${z(pH, s.p_hl2, s.n_trains).toFixed(2)})`);
check(minIsiOk, 'base neuron violated its own refractory period');
check(Math.abs(rateSum / s.n_trains - s.rate) < 0.01 * s.rate, 'realised rate');
check(Math.abs(contSum / s.n_trains - s.cont) < 0.005, 'realised contamination');
check(Math.abs(z(pS, s.p_sliding, s.n_trains)) < 3.5, 'Sliding RP acceptance matches Python');
check(Math.abs(z(pH, s.p_hl2, s.n_trains)) < 3.5, 'HL acceptance matches Python');

console.log(failures ? `\n${failures} FAILURE(S)` : '\nall checks passed');
process.exit(failures ? 1 : 0);
