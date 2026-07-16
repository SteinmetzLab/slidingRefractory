function [minCont, timeOfLowestCont] = computeMinContamination(...
    obsViol, spikeCount, refDur, rp, recDur, confThresh, rpReject, maxContam)
% COMPUTEMINCONTAMINATION  Analytical minimum confirmable contamination.
%
%   [minCont, timeOfLowestCont] = computeMinContamination(obsViol, spikeCount, ...
%       refDur, rp, recDur, confThresh, rpReject, maxContam)
%
%   Closed-form equivalent of scanning the contamination grid for the smallest
%   contamination that can be confirmed (rejected as the upper bound) at the
%   given confidence. For each tau_r the critical Poisson rate at which the
%   confidence equals confThresh is obtained from the Poisson<->chi-squared
%   identity:
%
%       lambda_crit = chi2inv(confThresh/100, 2*(obsViol + 1)) / 2
%
%   Inverting the Llobet expected-violations formula
%       Ve(C) = 2*tau/D * C*N * ((1-C)*N + (C*N - 1)/2)
%   for Ve = lambda_crit gives the smaller root
%
%       C_min(tau) = [ (N-0.5) - sqrt((N-0.5)^2 - lambda_crit*D/tau) ] / N.
%
%   The unit's minimum confirmable contamination is the minimum of C_min(tau)
%   over tau_r > rpReject. This matches the grid-based result of computeMatrix
%   to within the grid resolution, but is continuous and avoids building the
%   full confidence matrix.
%
%   INPUTS
%     obsViol    - [1 x nRP] cumulative observed violations V_o(tau_r).
%     spikeCount - total spike count N.
%     refDur     - [1 x nRP] tested tau_r (s) (right bin edge).
%     rp         - [1 x nRP] bin-centre tau_r (s), for the rpReject mask and
%                  the returned timeOfLowestCont.
%     recDur     - recording duration D (s).
%     confThresh - confidence threshold (%).
%     rpReject   - exclude tau_r at or below this (s).
%     maxContam  - (optional) cap in %; returns NaN above it. Default 35.
%
%   OUTPUTS
%     minCont          - minimum confirmable contamination (%), or NaN.
%     timeOfLowestCont - tau_r (s) at which it is achieved, or NaN.
%
%   SEE ALSO
%     slidingRP, computeViol, computeMatrix

if nargin < 8 || isempty(maxContam); maxContam = 35; end

N = spikeCount;
if N == 0 || recDur <= 0
    minCont = NaN; timeOfLowestCont = NaN; return;
end

gamma = confThresh / 100;
lambda = chi2inv(gamma, 2 * (obsViol + 1)) / 2;

disc = (N - 0.5)^2 - lambda .* recDur ./ refDur;
Cmin = nan(size(refDur));
ok = disc >= 0;
Cmin(ok) = ((N - 0.5) - sqrt(disc(ok))) / N * 100;   % percentage

testTimes = rp > rpReject;
CminTest = Cmin(testTimes);
if ~any(testTimes) || all(isnan(CminTest))
    minCont = NaN; timeOfLowestCont = NaN; return;
end

[minCont, idx] = min(CminTest);   % min ignores NaN by default
if minCont > maxContam
    minCont = NaN; timeOfLowestCont = NaN; return;
end
rpTest = rp(testTimes);
timeOfLowestCont = rpTest(idx);
end
