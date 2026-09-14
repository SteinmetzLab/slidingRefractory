function tests = test_power
% Unit tests for tauPass0 / minPassingFR and the appended slidingRP output.
%
%   results = runtests('matlab/tests/test_power.m')
%
% These quantities are diagnostics: they change no pass/fail decision. They
% answer "is there enough data for this test to say anything?", which the
% binary decision deliberately does not distinguish from "violations were
% observed". Reference values are cross-checked against the Python twins in
% python/slidingRP/power.py (tests/test_power.py).
tests = functiontests(localfunctions);
end

% =====================================================================
%  tauPass0
% =====================================================================

function test_tauPass0_reference_values(testCase)
% Values quoted in the manuscript revision (1 h recording, defaults:
% 10% contamination threshold, 90% confidence). Python agrees to 1e-12.
D = 3600;
verifyEqual(testCase, tauPass0(0.2*D, D)*1000, 84.22037648112823, 'RelTol', 1e-9);
verifyEqual(testCase, tauPass0(0.5*D, D)*1000, 13.4699, 'RelTol', 1e-4);
verifyEqual(testCase, tauPass0(1.1*D, D)*1000,  2.7823, 'RelTol', 1e-4);
end

function test_tauPass0_beyond_window_means_cannot_pass(testCase)
% A 0.2 spikes/s unit recorded for an hour needs a violation-free window of
% ~84 ms, far beyond the 10 ms the algorithm tests: it cannot pass whatever
% its autocorrelogram looks like.
D = 3600;
verifyGreaterThan(testCase, tauPass0(0.2*D, D), 0.010);
verifyLessThan(testCase,    tauPass0(2.0*D, D), 0.010);
end

function test_tauPass0_scaling(testCase)
% tauPass0 scales as D/N^2: doubling the rate quarters it, doubling the
% duration at fixed rate halves it.
D = 3600;
verifyEqual(testCase, tauPass0(1*D, D) / tauPass0(2*D, D), 4, 'RelTol', 0.02);
verifyEqual(testCase, tauPass0(1*D, D) / tauPass0(2*D, 2*D), 2, 'RelTol', 0.02);
end

function test_tauPass0_matches_grid_search(testCase)
% The closed form must agree with a direct search over the metric's own
% tau_r grid for the smallest tau at which zero violations passes.
D = 3600; acgBinSize = 1/30000;
for fr = [0.5 1 2 5]
    N = round(fr*D);
    rp = (0:acgBinSize:10/1000-acgBinSize) + acgBinSize/2;
    refDur = rp + acgBinSize/2;
    conf = 100 * computeViol(zeros(size(refDur)), [], N, refDur, 0.10, D);
    idx = find(conf >= 90, 1);
    if ~isempty(idx)
        % the grid lands on the first bin at or after the continuous answer
        d = refDur(idx) - tauPass0(N, D);
        verifyGreaterThanOrEqual(testCase, d, -1e-12);
        verifyLessThan(testCase, d, acgBinSize + 1e-12);
    else
        verifyGreaterThan(testCase, tauPass0(N, D), 10/1000);
    end
end
end

% =====================================================================
%  minPassingFR
% =====================================================================

function test_minPassingFR_inverts_tauPass0(testCase)
for D = [1800 3600 7200]
    for tau = [0.001 0.002 0.003 0.005]
        fr = minPassingFR(D, tau);
        verifyEqual(testCase, tauPass0(fr*D, D), tau, 'RelTol', 1e-9);
    end
end
end

function test_minPassingFR_matches_fig4g(testCase)
% Fig 4g: ~1.06 spikes/s at 1 h with a 3 ms RP, 90% confidence, 10% contamination.
verifyEqual(testCase, minPassingFR(3600, 0.003), 1.0593735747554274, 'RelTol', 1e-9);
% more confidence or a shorter clean window both demand more spikes
verifyGreaterThan(testCase, minPassingFR(3600, 0.003, 10, 95), minPassingFR(3600, 0.003));
verifyGreaterThan(testCase, minPassingFR(3600, 0.001), minPassingFR(3600, 0.003));
end

% =====================================================================
%  slidingRP appended output
% =====================================================================

function test_slidingRP_returns_tauPass0(testCase)
if isempty(which('histdiff'))
    assumeFail(testCase, 'histdiff (cortex-lab/spikes) not on path');
end
rng(1);
st = genST(5, 3600, 0.002);
[~, ~, ~, ~, ~, ~, ~, ~, ~, tp0] = slidingRP(st, struct('recDur', 3600));
verifyEqual(testCase, tp0, tauPass0(numel(st), 3600), 'RelTol', 1e-12);
end

function test_slidingRP_back_compatible_output_counts(testCase)
% Appending the output must not disturb existing call sites that request
% fewer outputs.
if isempty(which('histdiff'))
    assumeFail(testCase, 'histdiff (cortex-lab/spikes) not on path');
end
rng(2);
st = genST(5, 3600, 0.002);
p1 = slidingRP(st, struct('recDur', 3600));
[p4, c4] = slidingRP(st, struct('recDur', 3600));
[p9, c9, ~, ~, ~, ~, ~, ~, ~] = slidingRP(st, struct('recDur', 3600));
verifyEqual(testCase, p1, p4);
verifyEqual(testCase, p4, p9);
verifyEqual(testCase, c4, c9);
end

function test_slidingRP_tauPass0_in_corrected_path(testCase)
if isempty(which('histdiff'))
    assumeFail(testCase, 'histdiff (cortex-lab/spikes) not on path');
end
rng(3);
st = genST(5, 3600, 0.002);
st = st(1:min(5000, numel(st)));
params = struct('recDur', 3600, 'correction', true);
[~, ~, ~, ~, ~, ~, ~, ~, ~, tp0] = slidingRP(st, params);
verifyEqual(testCase, tp0, tauPass0(numel(st), 3600), 'RelTol', 1e-12);
end
