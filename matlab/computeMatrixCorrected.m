function[confMatrix, cont, rp, nACG, nominalConfMatrix] = computeMatrixCorrected(spikeTimes, params)

if nargin>1 && isfield(params, 'cont')
    cont = params.cont;
else
    cont = 0.5:0.5:35; 
end

if nargin>1 && isfield(params, 'recDur')
    recDur = params.recDur;
else
    recDur = max(spikeTimes); 
end

if nargin>1 && isfield(params, 'rpReject')
    rpReject = params.rpReject;
else
    rpReject = 0.0005;
end

rpBinSize = 1/30000;
rpEdges = 0:rpBinSize:10/1000;
spikeCount = numel(spikeTimes); 
firingRate = []; 

[nACG,rp] = histdiff(spikeTimes, spikeTimes, rpEdges);

confMatrix = nan(numel(cont), numel(rp));
nominalConfMatrix = nan(numel(cont), numel(rp));
expectedViolMatrix = nan(numel(cont), numel(rp));
obsViol = zeros(1, numel(rp));

validIdx = find(rp > rpReject);

% 1. Compute cumulative observed violations
for rpIdx = 1:numel(rp)
    obsViol(rpIdx) = sum(nACG(1,1:rpIdx));
end

% 2. Main computation loop over contamination levels
for cidx = 1:numel(cont)
    
    % A. Compute pointwise nominal confidences
    for rpIdx = 1:numel(rp)
        [nomScore, eViol] = computeViol(...
            obsViol(rpIdx), firingRate, spikeCount, rp(rpIdx)+rpBinSize/2, cont(cidx)/100, recDur);
        
        nominalConfMatrix(cidx, rpIdx) = nomScore * 100;
        expectedViolMatrix(cidx, rpIdx) = eViol;
    end
    
    % B. Multiple Comparisons Correction (Exact First-Passage Poisson DP)
    if isempty(validIdx)
        confMatrix(cidx, :) = 0;
        continue;
    end
    
    expectedViol = expectedViolMatrix(cidx, :);
    nomPvals = 1 - (nominalConfMatrix(cidx, :) / 100);
    minPval = min(nomPvals(validIdx));
    
    % Shortcut: Skip DP if it will fail anyway
    if minPval > 0.5
        confMatrix(cidx, :) = (1 - minPval) * 100;
        continue;
    end
    
    % --- OPTIMIZATION 1: Vectorized boundary finding ---
    c = -ones(1, numel(rp));
    expV_valid = expectedViol(validIdx);
    obsV_valid = obsViol(validIdx);
    max_obs = max(obsV_valid);
    
    if max_obs == 0
        % Fast-path for 0 observed violations
        c_valid = (exp(-expV_valid) <= minPval + 1e-10) - 1;
        c(validIdx) = c_valid;
    else
        % Explicitly match matrix sizes for poisscdf compatibility
        v_vec = (0:max_obs)';
        V_mat = repmat(v_vec, 1, numel(expV_valid));
        Exp_mat = repmat(expV_valid, max_obs + 1, 1);
        Obs_mat = repmat(obsV_valid, max_obs + 1, 1);
        
        % Vectorized computation of poisscdf for all v and all time bins
        cdf_matrix = poisscdf(V_mat, Exp_mat);
        
        % Logical mask of valid boundary conditions
        valid_mask = (V_mat <= Obs_mat) & (cdf_matrix <= minPval + 1e-10);
        
        % Because valid v's start at 0 and go up consecutively, the max valid v
        % is simply the total number of valid entries minus 1.
        c_valid = sum(valid_mask, 1) - 1;
        
        c(validIdx) = c_valid;
    end
    max_c = max(c(validIdx));
    
    % --- OPTIMIZATION 2: Vectorized Markov Chain DP ---
    if max_c < 0
        correctedConf = 1.0;
    else
        P_state = zeros(1, max_c + 1);
        P_state(1) = 1.0;
        P_false_pass = 0;
        
        % Pre-calculate lambda (expected additions) for all time steps
        lambda_all = [expectedViol(1), diff(expectedViol)];
        
        for k = 1:numel(rp)
            lambda_k = lambda_all(k);
            if lambda_k > 0
                pmf = poisspdf(0:max_c, lambda_k);
                
                % Use MATLAB's highly optimized built-in convolution
                P_next = conv(P_state, pmf);
                
                % Trim back to the tracked states
                P_state = P_next(1:max_c+1);
            end
            
            % Absorption: If it crossed the boundary c(k), it falsely passed
            if c(k) >= 0
                idx_cross = 1:(c(k) + 1);
                P_false_pass = P_false_pass + sum(P_state(idx_cross));
                P_state(idx_cross) = 0; % Remove from future tracking
            end
        end
        correctedConf = 1 - P_false_pass;
    end
    
    confMatrix(cidx, :) = correctedConf * 100;
end
end