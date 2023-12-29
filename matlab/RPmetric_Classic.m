

function [passTest, estContam, rp, nACG] = RPmetric_Classic(spikeTimes, params)

if nargin>1 && isfield(params, 'metricType')
    metricType = params.metricType;
else
    metricType = 'Llobet'; 
end
if nargin>1 && isfield(params, 'contaminationThresh')
    contaminationThresh = params.contaminationThresh;
else
    contaminationThresh = 10; % contamination level at which to test
end
if nargin>1 && isfield(params, 'recDur')
    recDur = params.recDur;
else
    recDur = max(spikeTimes); 
end
if nargin>1 && isfield(params, 'RPdur')
    RPdur = params.RPdur;
else
    RPdur = 0.002; % RP duration to test at 
end
if nargin>1 && isfield(params, 'nACG')
    % have provided a pre-calculated ACG. Must also have provided rp, the
    % bin labels of the ACG
    nACG = params.nACG;
    rp = params.rp;
    spikeCount = params.spikeCount;

    rpIdx = find(rp>RPdur,1);
    obsViol = sum(nACG(1:rpIdx));
    
else
    % calculate observed violations ourselves
    
    rpBinSize = 1/30000;
    rpEdges = 0:rpBinSize:RPdur;
    
    spikeCount = numel(spikeTimes); 
    
    [nACG,rp] = histdiff(spikeTimes, spikeTimes, rpEdges);
    
    obsViol = sum(nACG);

end

expectedNc = spikeCount * contaminationThresh/100;
expectedNb = spikeCount * (1-contaminationThresh/100); 

switch metricType

    % See paper for derivation of the expressions below. Note that the
    % "Hill" expressions are inaccurate but included for comparison. 

    case 'Llobet' % Llobet, Wyngaard, & Barbour 2022 - https://doi.org/10.1101/2022.02.08.479192
        %expectedViol = 2*P/D * N_c .* (N_b + (N_c-1)/2);
        expectedViol = 2*RPdur/recDur * expectedNc .* (expectedNb + (expectedNc-1)/2);
        %estContam = 1 - sqrt(1 - Nv*D/(Nt^2*P));
        estContam = 1 - sqrt(1 - obsViol*recDur/(spikeCount^2*RPdur));        

    case 'Hill' % Hill & Kleinfeld 2011 - https://doi.org/10.1523/JNEUROSCI.0971-11.2011
        %expectedViol = 2P/D * N_c * (N_c + N_b)
        expectedViol = 2*RPdur/recDur * expectedNc .* (expectedNb + expectedNc);
        %estContam = 1/2 * (1 - sqrt(1 - 2*Nv*D/Nt^2/P))
        estContam = 1/2 * (1 - sqrt(1 - 2*obsViol*recDur/spikeCount^2/RPdur));

end

passTest = obsViol<=expectedViol;
