function fr = minPassingFR(recDur, tau, varargin)
% MINPASSINGFR  Minimum firing rate for a violation-free unit to pass.
%
%   fr = minPassingFR(recDur, tau)
%   fr = minPassingFR(recDur, tau, contThresh, confThresh)
%
%   Exact inverse of Ve(tau) = -log(1 - gamma), which is a quadratic in the
%   spike count N:
%
%       [C(1-C) + C^2/2] N^2 - (C/2) N - log(1/(1-gamma)) * D / (2*tau) = 0
%
%   and returns N/D. This is the analytical form of manuscript Fig 4g, which
%   was previously obtained by simulation.
%
%   INPUTS
%     recDur     - recording duration in seconds.
%     tau        - refractory period duration assumed clean, in seconds.
%     contThresh - (optional) maximum acceptable contamination (%). Default 10.
%     confThresh - (optional) required confidence (%). Default 90.
%
%   OUTPUT
%     fr - minimum firing rate in spikes/s.
%
%   EXAMPLE
%     minPassingFR(3600, 0.003)   % 1.06 spikes/s at the paper defaults
%
%   SEE ALSO
%     tauPass0, slidingRP

if nargin > 2 && ~isempty(varargin{1}); contThresh = varargin{1}; else; contThresh = 10; end
if nargin > 3 && ~isempty(varargin{2}); confThresh = varargin{2}; else; confThresh = 90; end

C = contThresh / 100;
k = -log(1 - confThresh / 100) .* recDur ./ (2 * tau);
a = C * (1 - C) + C^2 / 2;
b = -C / 2;
N = (-b + sqrt(b.^2 + 4 * a * k)) / (2 * a);
fr = N ./ recDur;
end
