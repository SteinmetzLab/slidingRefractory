
function [passTest, estContam, rp, nACG] = RPmetric_Classic(spikeTimes, params)
% RPMETRIC_CLASSIC  Hill-Llobet point-estimate contamination metric (prior art).
%
%   [passTest, estContam, rp, nACG] = RPmetric_Classic(spikeTimes)
%   [passTest, estContam, rp, nACG] = RPmetric_Classic(spikeTimes, params)
%
%   Implements the prior state-of-the-art contamination metric as described
%   in Hill et al. (2011), with the Llobet et al. (2022) correction as the
%   default. Unlike the Sliding RP metric, this method uses a single fixed
%   RP duration and a point estimate of contamination (not a statistical
%   test). It is provided here for comparison with slidingRP().
%
%   The method has two known failure modes (see manuscript for details):
%     1. Always passes units with zero observed violations regardless of
%        firing rate (no statistical power check).
%     2. Incorrectly rejects clean units whose true RP is shorter than RPdur.
%
%   INPUTS
%     spikeTimes - [N x 1] vector of spike times in seconds.
%     params     - (optional) struct with any of the following fields:
%
%       .metricType         - 'Llobet' (default) or 'Hill'. Selects which
%                             expected-violations formula to use:
%                               'Llobet': Ve = 2P/D * Nc*(Nb + (Nc-1)/2)
%                               'Hill':   Ve = 2P/D * Nc*Nb
%       .contaminationThresh - Contamination threshold (%). Default: 10.
%       .RPdur              - Fixed assumed RP duration (seconds).
%                             Default: 0.002 (2 ms).
%       .recDur             - Recording duration (seconds). Defaults to
%                             max(spikeTimes) if not specified.
%       .acgBinSize         - ACG bin width (seconds). Default: 1/30000.
%       .nACG               - Pre-computed ACG counts (skips histdiff).
%                             Must also provide .rp and .spikeCount.
%       .rp                 - Bin labels for pre-computed nACG (seconds).
%       .spikeCount         - Total spike count for pre-computed ACG.
%
%   OUTPUTS
%     passTest  - Logical. True if obsViol <= expectedViol at RPdur.
%     estContam - Point estimate of contamination proportion (not %).
%                 Derived by inverting the expected-violations formula.
%     rp        - Bin label vector from ACG computation (seconds).
%     nACG      - ACG counts used for the computation.
%
%   REFERENCES
%     Hill DN, Mehta SB, Kleinfeld D (2011) Quality metrics to accompany
%       spike sorting of extracellular signals. J Neurosci 31:8699-8705.
%     Llobet V, Wyngaard A, Barbour B (2022) Automatic post-processing and
%       merging of multiple spike-sorting analyses. bioRxiv.

if nargin > 1 && isfield(params, 'metricType')
    metricType = params.metricType;
else
    metricType = 'Llobet';
end
if nargin > 1 && isfield(params, 'contaminationThresh')
    contaminationThresh = params.contaminationThresh;
else
    contaminationThresh = 10;
end
if nargin > 1 && isfield(params, 'recDur')
    recDur = params.recDur;
else
    recDur = max(spikeTimes);
end
if nargin > 1 && isfield(params, 'RPdur')
    RPdur = params.RPdur;
else
    RPdur = 0.002;
end
if nargin > 1 && isfield(params, 'acgBinSize')
    acgBinSize = params.acgBinSize;
else
    acgBinSize = 1/30000;
end

if nargin > 1 && isfield(params, 'nACG')
    % Use pre-computed ACG (avoids histdiff dependency)
    nACG = params.nACG;
    rp   = params.rp;
    spikeCount = params.spikeCount;
    rpIdx = find(rp > RPdur, 1);
    obsViol = sum(nACG(1:rpIdx));
else
    rpEdges    = 0:acgBinSize:RPdur;
    spikeCount = numel(spikeTimes);
    [nACG, rp] = histdiff(spikeTimes, spikeTimes, rpEdges);
    obsViol    = sum(nACG);
end

expectedNc = spikeCount * contaminationThresh / 100;
expectedNb = spikeCount * (1 - contaminationThresh / 100);

switch metricType

    case 'Llobet'
        % Llobet, Wyngaard & Barbour (2022): contaminating spikes also
        % produce violations with each other.
        %   Ve = 2P/D * Nc * (Nb + (Nc-1)/2)
        %   C  = 1 - sqrt(1 - Nv*D / (Nt^2 * P))
        expectedViol = 2 * RPdur / recDur * expectedNc .* (expectedNb + (expectedNc - 1) / 2);
        estContam    = 1 - sqrt(1 - obsViol * recDur / (spikeCount^2 * RPdur));

    case 'Hill'
        % Hill, Mehta & Kleinfeld (2011): contamination from a single other
        % neuron, which produces violations only with the base neuron (not with
        % itself), so Ve uses Nb (not Nb + Nc). Matches the manuscript Methods.
        %   Ve = 2P/D * Nc * Nb
        %   C  = 1/2 * (1 - sqrt(1 - 2*Nv*D / (Nt^2 * P)))
        expectedViol = 2 * RPdur / recDur * expectedNc .* expectedNb;
        estContam    = 1/2 * (1 - sqrt(1 - 2 * obsViol * recDur / spikeCount^2 / RPdur));

end

passTest = obsViol <= expectedViol;
