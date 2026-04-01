
function plotSlidingRP(spikeTimes, params)
% PLOTSLIDINGRP  Visualise the Sliding RP metric result for one cluster.
%
%   plotSlidingRP(spikeTimes)
%   plotSlidingRP(spikeTimes, params)
%
%   Produces a three-panel figure:
%     Panel 1 – Autocorrelogram (ACG) from 0 to 5 ms.
%     Panel 2 – Full confidence matrix (contamination vs tau_r), with the
%               90%-confidence iso-contour and the identified minimum-
%               contamination time overlaid.
%     Panel 3 – Confidence trace at the 10% contamination threshold.
%
%   Title text is coloured green (pass), red (fail with violations), or
%   blue (fail due to low firing rate, i.e. zero short-ISI spikes).
%
%   INPUTS
%     spikeTimes - [N x 1] vector of spike times in seconds.
%     params     - (optional) struct passed through to slidingRP(). Any
%                  field accepted by slidingRP() may be used. Additional
%                  field:
%                    .cidx - Cluster ID to display in the figure title.
%
%   SEE ALSO
%     slidingRP, slidingRP_all

if nargin < 2 || ~isstruct(params)
    params = struct(); params.cidx = [];
else
    if ~isfield(params, 'cidx'); params.cidx = []; end
end

[passTest, maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont, ...
    nSpikesBelow2, confMatrix, cont, rp, nACG] = slidingRP(spikeTimes, params);

if isfield(params, 'recDur'); recDur = params.recDur; else; recDur = max(spikeTimes); end
firingRate = numel(spikeTimes) / recDur;

f = figure; f.Color = 'w';
fp = f.Position;
f.Position = [fp(1) fp(2) 1300 369];

% --- Panel 1: ACG ---
subplot(1, 3, 1);
bar(rp * 1000, nACG, 1, 'FaceColor', 'k', 'EdgeAlpha', 0);
xlim([0 5]);
xlabel('Time from spike (ms)');
ylabel('ACG count (spks)');
if ~isempty(params.cidx)
    t1 = title(sprintf('Cluster #%d: FR=%.2f', params.cidx, firingRate));
else
    t1 = title(sprintf('FR=%.2f', firingRate));
end
hold on;
fill([0 1 1 0] * 0.5, [0 0 1 1] * max(ylim()), 'k', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
box off;

% --- Panel 2: Confidence matrix ---
subplot(1, 3, 2); hold on;
imagesc(rp * 1000, cont, confMatrix);
plot([rp(1) rp(end)] * 1000, [10 10], 'r');   % contamination threshold line

if ~isnan(timeOfLowestCont)
    plot(timeOfLowestCont * 1000 * [1 1], [cont(1) cont(end)], 'r');

    % 90%-confidence iso-contour
    [~, ii] = max([zeros(1, numel(rp)); confMatrix] > 90, [], 1);
    ii(ii == 1) = nan;
    contContour = nan(size(ii));
    contContour(~isnan(ii)) = cont(ii(~isnan(ii)) - 1);
    plot(rp * 1000, contContour, 'r', 'LineWidth', 2.0);
end

fill([0 1 1 0] * 0.5, [0 0 1 1] * max(ylim()), 'k', 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
h = colorbar;
h.Label.String = 'Confidence (%)';
caxis([0 100]);
xlabel('Time from spike (ms)');
xlim([0 5]);
ylabel('Contamination (%)');
ylim([0 max(cont)]);
set(gca, 'YDir', 'reverse');
t2 = title(sprintf('max conf = %.2f%%, min cont = %.1f%%, time = %.2f ms', ...
    maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont * 1000));

% --- Panel 3: Confidence at 10% contamination ---
subplot(1, 3, 3);
plot(rp * 1000, confMatrix(cont == 10, :), 'k', 'LineWidth', 2.0);
xlabel('Time from spike (ms)');
ylabel('Confidence of \leq10% contamination (%)');
box off; hold on;
plot([0 5], [90 90], 'r');
fill([0 1 1 0] * 0.5, [0 0 1 1] * 100, 'k', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
xlim([0 5]);
ylim([0 100]);

% Colour titles: green=pass, blue=fail (no violations), red=fail
if minContWith90Confidence <= 10
    t1.Color = [34 177 76] / 255;
    t2.Color = [34 177 76] / 255;
elseif nSpikesBelow2 == 0
    t1.Color = 'b';
    t2.Color = 'b';
else
    t1.Color = 'r';
    t2.Color = 'r';
end
