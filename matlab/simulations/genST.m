

function st = genST(rate, duration, varargin)
% GENST  Generate a simulated spike train (Poisson, optionally with a hard RP).
%
%   st = genST(rate, duration)
%   st = genST(rate, duration, refractoryPeriod)
%
%   Draws inter-spike intervals as refractoryPeriod + Exponential(mu) and
%   returns their cumulative sum (spike times). The exponential mean mu is
%   rate-corrected so that, despite the dead time imposed by the refractory
%   period, the realised total rate still equals `rate`:
%
%       rSim = rate / (1 - refractoryPeriod*rate);   mu = 1/rSim
%
%   This matches the simulation in Roth et al. (2026), Methods ("Spike train
%   simulations"): I = P + exprnd(mu), mu = (1 - P*R)/R. To build a
%   contaminated train, call once for the base neuron (with its RP) and once
%   for the contamination (refractoryPeriod = 0), then concatenate and sort.
%
%   INPUTS
%     rate            - Desired total firing rate (same time units as duration).
%     duration        - Recording duration.
%     refractoryPeriod - (optional) Hard refractory period; default 0.
%
%   OUTPUT
%     st - Column vector of spike times in [0, duration), in the same units.
%
%   NOTE: spike times are continuous reals, not quantised to a sample grid.
%   See also genST.m in python/slidingRP/simulations.py for the Python twin.

if nargin<3
    refPeriod = 0;
else
    refPeriod = varargin{1};
end

%%

rSim = rate/(1-refPeriod*rate); % this adjustment makes it so that the spikes generated with a RP will have the correct total rate

mu = 1/rSim; % mean of distribution of ISIs

n = rate*duration; % number of spikes we'll need

isi = refPeriod+exprnd(mu, ceil(n*2), 1); % generating twice as many spikes as we expect, we'll cut the extras

while sum(isi)<duration
    % this part will be extremely slow if it needs to be invoked, but in
    % general it should not be - only if we didn't get enough spikes
    % somehow
    isi(end+1) = exprnd(mu);
end

st = cumsum(isi); % convert to spike times
st = st(1:find(st<duration,1,'last')); %cut all the ones that we don't need

return;

% ------------------------------------------------------------------------
% The cells below are unreachable (after the return above). They are kept as
% manual sanity checks: copy a cell into the command window to confirm the
% realised rate matches the requested rate and the ISI/ACG look correct.
% ------------------------------------------------------------------------

%% test code

rate = 5; n = 10000;
st = genST(rate,n);

sc = hist(st,0:n); % count spikes in 1 sec intervals
figure; 
x = 0:20;
hist(sc,x);
hold on; 
plot(x, poisspdf(x,rate)*n, 'r', 'LineWidth', 2.0);

%% test 2

rate = 20; n = 100000; rp = 0.004; 
st = genST(rate,n, rp); 

sc = hist(st,0:n); % count spikes in 1 sec intervals
figure; 

subplot(1,2,1); 
x = 0:max(sc)+2;
hist(sc,x);
hold on; 
plot(x, poisspdf(x,rate)*n, 'ro-', 'LineWidth', 2.0);
title(sprintf('mean = %.3f', mean(sc)));

subplot(1,2,2); 
 [xLin, nLin, xLog, nLog] = myACG(st, gca, []);
xlim([0 0.02]);
