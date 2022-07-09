

function plotSlidingRP(spikeTimes, params)

[maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,...
    nSpikesBelow2, confMatrix, cont, rp, nACG, firingRate] ...
    = slidingRP(spikeTimes, params);

f = figure; f.Color = 'w';
fp = f.Position;
f.Position = [fp(1) fp(2) 1300 369];

subplot(1,3,1); 
bar(rp*1000, nACG, 1,'FaceColor', 'k', 'EdgeAlpha', 0); 
xlim([0 5]); 
xlabel('Time from spike (ms)');
ylabel('ACG count (spks)'); 
t1 = title(sprintf('Cluster #%d: FR=%.2f', params.cidx, firingRate)); 
hold on; 
fill([0 1 1 0]*0.5, [0 0 1 1]*max(ylim()), 'k', 'FaceAlpha', 0.2, 'EdgeAlpha', 0); 
box off; 


subplot(1,3,2); hold on;
imagesc(rp*1000, cont, confMatrix);
plot([rp(1) rp(end)]*1000, [10 10], 'r'); 

if ~isnan(timeOfLowestCont)
    plot(timeOfLowestCont*1000*[1 1], [cont(1) cont(end)],'r');

    % compute the conf=90 contour
    [~,ii] = max([zeros(1, numel(rp)); confMatrix]>90, [], 1);
    ii(ii==1) = nan;
    contContour = nan(size(ii)); contContour(~isnan(ii)) = cont(ii(~isnan(ii))-1);
    plot(rp*1000, contContour, 'r', 'LineWidth', 2.0); 
end

fill([0 1 1 0]*0.5, [0 0 1 1]*max(ylim()), 'k', 'FaceAlpha', 0.5, 'EdgeAlpha', 0); 

h = colorbar;
h.Label.String = 'Confidence (%)'; 
caxis([0 100]); 
xlabel('Time from spike (ms)');
xlim([0 5]); 
ylabel('Contamination (%)'); 
ylim([0 max(cont)]); 
set(gca, 'YDir', 'reverse'); 
t2 = title(sprintf('max conf = %.2f%%, min cont = %.1f%%, time = %.2f ms', ...
    maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont*1000));


if minContWith90Confidence <= 10
    t1.Color = [34 177 76]/255;
    t2.Color = [34 177 76]/255;
elseif nSpikesBelow2==0
    t1.Color = 'b';
    t2.Color = 'b';
else
    t1.Color = 'r';
    t2.Color = 'r';
end

subplot(1,3,3); 
plot(rp*1000, confMatrix(cont==10,:), 'k','LineWidth', 2.0)
xlabel('Time from spike (ms)');
ylabel('Confidence of \leq10% contamination (%)'); 
box off; hold on; 
plot([0 5], [90 90], 'r'); 
fill([0 1 1 0]*0.5, [0 0 1 1]*100, 'k', 'FaceAlpha', 0.2, 'EdgeAlpha', 0); 
xlim([0 5]); 
ylim([0 100]); 