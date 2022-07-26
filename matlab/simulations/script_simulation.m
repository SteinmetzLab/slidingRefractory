

%% simulation test


clear all; 

nSim = 100; 

binSize = 1/30000;
b = (1:300)/30000;

baseRates = logspace(-0.3,1.3,8); 
recDur = 3600; 
rp = 0.002;
contPct = 0:0.02:0.2;
params = struct(); params.cont = 10;

figure; 
passPct = zeros(numel(baseRates), numel(contPct));
colors = hsv(numel(baseRates)); 
for bidx = 1:numel(baseRates)
    baseRate = baseRates(bidx);
    for c = 1:numel(contPct)
        
        contRate = baseRate*contPct(c);
        
        
        simRes = zeros(nSim, 1);
        for n = 1:nSim
            
            st = genST(baseRate,recDur);
            isi = diff([0; st]); isi(isi<rp) = [];
            st = cumsum(isi);
            
            contST = genST(contRate,recDur);
            
            combST = sort([st; contST]);
 
            [confMatrix, cont, rpTestVals, ~, ~] = computeMatrix(combST, params);
            simRes(n) = max(confMatrix(rpTestVals>0.0005))>90;
        end
        
        passPct(bidx,c) = sum(simRes)/nSim*100;
        
        
    end
    legH(bidx) = plot(contPct, passPct(bidx,:), 'o-', 'Color', colors(bidx,:), 'MarkerFaceColor', colors(bidx,:)); 
        
    hold on; drawnow;
end
addX(0.1);
legend(legH, array2stringCell(baseRates));
set(gcf, 'Color', 'w'); 
box off; 
xlabel('True contamination proportion'); 
ylabel('Percentage of simulations that pass'); 


%% test the 'expected count'

baseRates = logspace(-0.6,1.5,80); 
recDur = 7200; 
rp = 0.002;
contProp = 0.1;
params = struct(); params.cont = 10;

for bidx = 1:numel(baseRates)
    baseRate = baseRates(bidx);
    contRate = baseRate*contProp;
    
    st = genST(baseRate,recDur);
    isi = diff([0; st]); isi(isi<rp) = [];
    st = cumsum(isi);
    
    contST = genST(contRate,recDur);
    
    combST = sort([st; contST]);
    
    [confMatrix, cont, rpTestVals, nACG, firingRate] = computeMatrix(combST, params);
    
    contaminationRate = firingRate*contProp; 
    expectedViol = contaminationRate*rp*2*numel(combST);
    
    obsExp(bidx,1) = sum(nACG(rpTestVals<=0.002)); 
    obsExp(bidx,2) = expectedViol;
end

f = figure; f.Color = 'w';
plot(baseRates, obsExp(:,1), 'r.-'); hold on;
plot(baseRates, obsExp(:,2), 'k.-');
xlabel('firing rate (sp/s)')
ylabel('count'); 
legend({'observed', 'expected'});