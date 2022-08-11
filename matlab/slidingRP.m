

function [maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,...
    nSpikesBelow2, confMatrix, cont, rp, nACG, firingRate] ...
    = slidingRP(spikeTimes, params)
% compute the metric for a single cluster in a recording

% returns:
% - Max confidence that you have <= 10% contamination
% - Minimum contamination for which you have >=90% confidence
% - Time at which best score happens



[confMatrix, cont, rp, nACG, firingRate] = computeMatrix(spikeTimes, params); 
% matrix is [nCont x nRP]

testTimes = rp>0.0005; 

maxConfidenceAt10Cont = max(confMatrix(cont==10, testTimes)); 

[ii,jj] = find(confMatrix(:,testTimes)>90); 
[minI, idx] = min(ii); 
minContWith90Confidence = cont(minI);  
if isempty(minContWith90Confidence); minContWith90Confidence = NaN; end

[~,minRP] = max(confMatrix(minI,testTimes)); 
timeOfLowestCont = rp(minRP+find(testTimes,1));
if isempty(timeOfLowestCont); timeOfLowestCont = NaN; end

nSpikesBelow2 = sum(nACG(1:find(rp>0.002,1)));