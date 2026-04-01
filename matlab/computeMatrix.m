
function [confMatrix, cont, rp, nACG] = computeMatrix(spikeTimes, params)
% COMPUTEMATRIX  Build the confidence matrix for the Sliding RP metric.
%
%   [confMatrix, cont, rp, nACG] = computeMatrix(spikeTimes)
%   [confMatrix, cont, rp, nACG] = computeMatrix(spikeTimes, params)
%
%   For each combination of candidate RP duration (tau_r) and contamination
%   level (C), computes the statistical confidence that the observed number
%   of ACG violations is consistent with contamination < C. The result is a
%   2-D confidence matrix of size [nCont x nRP].
%
%   Requires histdiff from cortex-lab/spikes on the MATLAB path.
%
%   INPUTS
%     spikeTimes - [N x 1] vector of spike times in seconds.
%     params     - (optional) struct with any of the following fields:
%
%       .recDur      - Recording duration (seconds). Defaults to
%                      max(spikeTimes) if not specified. Recommended to
%                      set explicitly for accurate expected-violation counts.
%       .cont        - Vector of contamination levels (%) to test.
%                      Default: 0.5:0.5:35 (70 levels).
%       .acgBinSize  - ACG bin width (seconds). Default: 1/30000 s
%                      (one sample at 30 kHz).
%       .testWindow  - Maximum tau_r to test (seconds). Default: 0.01 (10 ms).
%
%   OUTPUTS
%     confMatrix - [nCont x nRP] matrix. Entry (i, j) is the confidence (%)
%                  that contamination is below cont(i) when testing at
%                  tau_r = rp(j). Computed via computeViol().
%     cont       - [1 x nCont] vector of contamination levels tested (%).
%     rp         - [1 x nRP] vector of tau_r bin centres tested (seconds).
%                  Each value is the left edge of the ACG bin as returned
%                  by histdiff; the right edge used in computeViol is
%                  rp(j) + acgBinSize/2.
%     nACG       - [1 x nRP] ACG counts: number of spike pairs with ISI
%                  in each bin.
%
%   SEE ALSO
%     slidingRP, computeViol

if nargin > 1 && isfield(params, 'cont')
    cont = params.cont;
else
    cont = 0.5:0.5:35;
end
if nargin > 1 && isfield(params, 'recDur')
    recDur = params.recDur;
else
    recDur = max(spikeTimes);
end
if nargin > 1 && isfield(params, 'acgBinSize')
    acgBinSize = params.acgBinSize;
else
    acgBinSize = 1/30000;
end
if nargin > 1 && isfield(params, 'testWindow')
    testWindow = params.testWindow;
else
    testWindow = 0.01;
end

rpEdges = 0:acgBinSize:testWindow;

spikeCount = numel(spikeTimes);
firingRate = []; % unused in current Llobet formulation; kept for API compatibility

[nACG, rp] = histdiff(spikeTimes, spikeTimes, rpEdges);
% nACG(1, k): number of spike pairs with ISI in bin k
% rp(k): left edge of bin k (seconds)

confMatrix = nan(numel(cont), numel(rp));
for rpIdx = 1:numel(rp)

    % Cumulative observed violations up to tau_r = rp(rpIdx)
    obsViol = sum(nACG(1, 1:rpIdx));

    for cidx = 1:numel(cont)
        % Use right edge of bin (rp + half bin width) as the tested tau_r
        confMatrix(cidx, rpIdx) = 100 * computeViol(...
            obsViol, firingRate, spikeCount, ...
            rp(rpIdx) + acgBinSize/2, cont(cidx)/100, recDur);
    end
end
