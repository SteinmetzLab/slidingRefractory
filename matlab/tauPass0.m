function tau = tauPass0(spikeCount, recDur, varargin)
% TAUPASS0  Shortest violation-free window that would let a unit pass.
%
%   tau = tauPass0(spikeCount, recDur)
%   tau = tauPass0(spikeCount, recDur, contThresh, confThresh)
%
%   A unit can fail the Sliding RP test for two different reasons: refractory
%   period violations incompatible with acceptable contamination, or too few
%   spikes to establish acceptable contamination at all. The pass/fail decision
%   does not distinguish them (what counts as a plausible refractory window is
%   the user's assumption, not the method's), but this quantity makes the
%   distinction explicit.
%
%   Under the Llobet model the expected violation count up to tau_r is
%       Ve = 2*tau_r*Nc*(Nb + (Nc-1)/2)/D,   Nc = C*N,  Nb = (1-C)*N
%   so a unit with ZERO observed violations up to tau_r passes at confidence
%   gamma iff 1 - exp(-Ve) >= gamma. Inverting for tau_r:
%
%       tauPass0 = -log(1 - gamma) * D / (2*Nc*(Nb + (Nc-1)/2))
%
%   A unit whose tauPass0 exceeds the longest tested RP duration (10 ms by
%   default) cannot pass however clean its autocorrelogram is: there is not
%   enough data. One whose tauPass0 is below a refractory duration the user
%   considers plausible has enough data for the test to be informative.
%
%   INPUTS
%     spikeCount - total number of spikes in the unit (scalar or array).
%     recDur     - recording duration in seconds (scalar or array).
%     contThresh - (optional) maximum acceptable contamination (%). Default 10.
%     confThresh - (optional) required confidence (%). Default 90.
%
%   OUTPUT
%     tau - shortest violation-free window in seconds; Inf if spikeCount is 0.
%
%   EXAMPLE (defaults, 1 hour recording)
%     tauPass0(0.2*3600, 3600)*1000   % 84.2 ms -- cannot pass at all
%     tauPass0(1.1*3600, 3600)*1000   %  2.8 ms -- evaluable for a 3 ms RP
%
%   SEE ALSO
%     slidingRP, minPassingFR

if nargin > 2 && ~isempty(varargin{1}); contThresh = varargin{1}; else; contThresh = 10; end
if nargin > 3 && ~isempty(varargin{2}); confThresh = varargin{2}; else; confThresh = 90; end

C  = contThresh / 100;
Nc = spikeCount .* C;
Nb = spikeCount .* (1 - C);
denom = 2 * Nc .* (Nb + (Nc - 1) / 2);

tau = -log(1 - confThresh / 100) .* recDur ./ denom;
tau(denom <= 0) = Inf;
end
