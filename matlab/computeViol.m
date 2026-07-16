
function [confidenceScore, expectedViol] = computeViol(...
    obsViol, ~, spikeCount, refDur, contaminationProp, recDur)
% Second argument (firingRate) is accepted but unused; pass [] for clarity.
% COMPUTEVIOL  Compute the Poisson confidence score for a single (tau_r, C) pair.
%
%   [confidenceScore, expectedViol] = computeViol(obsViol, firingRate, ...
%       spikeCount, refDur, contaminationProp, recDur)
%
%   For a hypothesised contamination proportion C, computes the expected
%   number of refractory period (RP) violations under the Llobet et al.
%   (2022) model and returns the probability that the observed violation
%   count arose from a unit with contamination >= C (i.e. the confidence
%   that contamination is *less than* C).
%
%   The expected violations formula (Llobet et al. 2022) accounts for
%   contaminating spikes producing violations both with base-neuron spikes
%   and with each other:
%
%     Ve = (2 * tau_r / D) * Nc * (Nb + (Nc - 1) / 2)
%
%   where Nc = C * N_total  (contaminating spikes)
%         Nb = (1-C) * N_total  (base-neuron spikes)
%         D  = recording duration
%
%   The confidence score is then:
%
%     gamma = 1 - Poisson_CDF(obsViol ; Ve)
%           = P(V > obsViol | lambda = Ve)
%
%   Intuitively: if we observe fewer violations than expected for
%   contamination level C, gamma is high (confident the unit is clean).
%
%   INPUTS
%     obsViol          - Observed number of RP violations (cumulative ACG
%                        count up to tau_r). Scalar or vector.
%     firingRate       - (unused, kept for API compatibility) Pass [].
%     spikeCount       - Total number of spikes in the spike train (N_total).
%     refDur           - Refractory period duration tested, tau_r (seconds).
%                        Typically rp(k) + acgBinSize/2, i.e. the right edge
%                        of ACG bin k.
%     contaminationProp - Hypothesised contamination as a proportion in [0,1]
%                        (e.g. 0.10 for 10%). NOT a percentage.
%     recDur           - Recording duration D (seconds).
%
%   OUTPUTS
%     confidenceScore  - Probability in [0, 1] that true contamination is
%                        below contaminationProp given the observed violation
%                        count. Multiply by 100 for percentage.
%     expectedViol     - Expected number of violations Ve under the
%                        hypothesised contamination level.
%
%   REFERENCE
%     Llobet V, Wyngaard A, Barbour B (2022) Automatic post-processing and
%     merging of multiple spike-sorting analyses. bioRxiv.

% Llobet et al. 2022 formulation
Nc = spikeCount * contaminationProp;
Nb = spikeCount * (1 - contaminationProp);
expectedViol = 2 * refDur / recDur * Nc .* (Nb + (Nc - 1) / 2);

confidenceScore = 1 - poisscdf(obsViol, expectedViol);
