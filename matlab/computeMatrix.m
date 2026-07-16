
function [confMatrix, cont, rp, nACG, nominalConfMatrix] = computeMatrix(spikeTimes, params)
% COMPUTEMATRIX  Build the confidence matrix for the Sliding RP metric.
%
%   [confMatrix, cont, rp, nACG] = computeMatrix(spikeTimes)
%   [confMatrix, cont, rp, nACG] = computeMatrix(spikeTimes, params)
%   [confMatrix, cont, rp, nACG, nominalConfMatrix] = computeMatrix(...)
%
%   For each combination of candidate RP duration (tau_r) and contamination
%   level (C), computes the statistical confidence that the observed number
%   of ACG violations is consistent with contamination < C. The result is a
%   2-D confidence matrix of size [nCont x nRP].
%
%   By default the confidence is the pointwise (per-tau_r) value. Setting
%   params.correction = true instead applies a family-wise-error-rate
%   correction across the swept tau_r values (an exact Poisson first-passage /
%   Markov-chain calculation; see Methods, "False acceptance rate in the
%   Sliding RP metric", Fig. S3). The correction is computationally expensive
%   and, because it cannot be exact without knowing the true RP duration,
%   over-penalises short-RP units; it is provided for Fig. S3 only and is off
%   by default.
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
%       .correction  - Logical. Apply the multiple-comparisons correction.
%                      Default: false.
%       .rpReject    - Min tau_r (s) included in the correction. Default
%                      0.0005. (Only used when correction is true.)
%
%   OUTPUTS
%     confMatrix - [nCont x nRP] matrix. Entry (i, j) is the confidence (%)
%                  that contamination is below cont(i) when testing at
%                  tau_r = rp(j). Pointwise, or FWER-corrected when
%                  params.correction is true.
%     cont       - [1 x nCont] vector of contamination levels tested (%).
%     rp         - [1 x nRP] vector of tau_r bin centres tested (seconds).
%                  Each value is the left edge of the ACG bin as returned
%                  by histdiff; the right edge used in computeViol is
%                  rp(j) + acgBinSize/2.
%     nACG       - [1 x nRP] ACG counts: number of spike pairs with ISI
%                  in each bin.
%     nominalConfMatrix - [nCont x nRP] uncorrected (pointwise) confidence (%).
%                  Equal to confMatrix when correction is false; returned for
%                  comparison when correction is true.
%
%   SEE ALSO
%     slidingRP, computeViol, computeMinContamination

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
if nargin > 1 && isfield(params, 'correction')
    correction = params.correction;
else
    correction = false;
end
if nargin > 1 && isfield(params, 'rpReject')
    rpReject = params.rpReject;
else
    rpReject = 0.0005;
end

rpEdges = 0:acgBinSize:testWindow;

spikeCount = numel(spikeTimes);

[nACG, rp] = histdiff(spikeTimes, spikeTimes, rpEdges);
% nACG(1, k): number of spike pairs with ISI in bin k
% rp(k): left edge of bin k (seconds)

obsViol = cumsum(nACG);             % cumulative observed violations up to tau_r
refDur  = rp + acgBinSize/2;        % right edge of bin = tested tau_r

% ---- Pointwise (nominal) confidence matrix, computed for every case --------
nominalConfMatrix    = nan(numel(cont), numel(rp));
expectedViolMatrix   = nan(numel(cont), numel(rp));
for cidx = 1:numel(cont)
    [nomScore, eViol] = computeViol(...
        obsViol, [], spikeCount, refDur, cont(cidx)/100, recDur);
    nominalConfMatrix(cidx, :)  = nomScore * 100;
    expectedViolMatrix(cidx, :) = eViol;
end

if ~correction
    confMatrix = nominalConfMatrix;
    return;
end

% ---- Multiple-comparisons correction (exact Poisson first-passage DP) -------
% For each contamination level, compute the probability that a random Poisson
% accumulation crosses the per-bin passing boundary anywhere across the swept
% tau_r window; the corrected confidence is 1 - that false-passing probability.
confMatrix = nan(numel(cont), numel(rp));
validIdx = find(rp > rpReject);

for cidx = 1:numel(cont)

    if isempty(validIdx)
        confMatrix(cidx, :) = 0;
        continue;
    end

    expectedViol = expectedViolMatrix(cidx, :);
    nomPvals = 1 - (nominalConfMatrix(cidx, :) / 100);
    minPval = min(nomPvals(validIdx));

    % Shortcut: if even the best pointwise p-value is poor, the corrected
    % confidence cannot exceed (1 - minPval); skip the DP.
    if minPval > 0.5
        confMatrix(cidx, :) = (1 - minPval) * 100;
        continue;
    end

    % --- Vectorised boundary finding: max observed count per bin that still
    %     yields a pointwise p-value <= minPval ---
    c = -ones(1, numel(rp));
    expV_valid = expectedViol(validIdx);
    obsV_valid = obsViol(validIdx);
    max_obs = max(obsV_valid);

    if max_obs == 0
        c_valid = (exp(-expV_valid) <= minPval + 1e-10) - 1;
        c(validIdx) = c_valid;
    else
        v_vec   = (0:max_obs)';
        V_mat   = repmat(v_vec, 1, numel(expV_valid));
        Exp_mat = repmat(expV_valid, max_obs + 1, 1);
        Obs_mat = repmat(obsV_valid, max_obs + 1, 1);
        cdf_matrix = poisscdf(V_mat, Exp_mat);
        valid_mask = (V_mat <= Obs_mat) & (cdf_matrix <= minPval + 1e-10);
        c_valid = sum(valid_mask, 1) - 1;
        c(validIdx) = c_valid;
    end
    max_c = max(c(validIdx));

    % --- Markov-chain DP for the first-passage probability ---
    if max_c < 0
        correctedConf = 1.0;
    else
        P_state = zeros(1, max_c + 1);
        P_state(1) = 1.0;
        P_false_pass = 0;
        lambda_all = [expectedViol(1), diff(expectedViol)];

        for k = 1:numel(rp)
            lambda_k = lambda_all(k);
            if lambda_k > 0
                pmf = poisspdf(0:max_c, lambda_k);
                P_next = conv(P_state, pmf);
                P_state = P_next(1:max_c + 1);
            end
            if c(k) >= 0
                idx_cross = 1:(c(k) + 1);
                P_false_pass = P_false_pass + sum(P_state(idx_cross));
                P_state(idx_cross) = 0;
            end
        end
        correctedConf = 1 - P_false_pass;
    end

    confMatrix(cidx, :) = correctedConf * 100;
end
end
