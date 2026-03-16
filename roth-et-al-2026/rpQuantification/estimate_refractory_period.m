function [RP_estimated, sigmoid_params, fit_quality] = estimate_refractory_period(acg, bin_centers, bin_size, P_slope)
% Estimate refractory period from single-sided autocorrelogram using sigmoid fit
%
% Inputs:
%   acg - autocorrelogram values (firing rates), single-sided starting at t=0
%   bin_centers - time values for each bin (in seconds), should start at or near 0
%   bin_size - size of each bin (in seconds), default: 1/30000
%   P_slope - percent of slope for RP estimation, default: 0.1
%
% Outputs:
%   RP_estimated - estimated refractory period in milliseconds
%   sigmoid_params - parameters of fitted sigmoid [ymin, ymax, k, x0]
%   fit_quality - struct with goodness of fit metrics

    % Set defaults
    if nargin < 3 || isempty(bin_size)
        bin_size = 1/30000;
    end
    if nargin < 4 || isempty(P_slope)
        P_slope = 0.1;
    end
    
    % Validate input
    if min(bin_centers) < 0
        warning('ACG should be single-sided (starting at t>=0). Using absolute values.');
        bin_centers = abs(bin_centers);
        [bin_centers, sort_idx] = sort(bin_centers);
        acg = acg(sort_idx);
    end
    
    % Initialize outputs
    RP_estimated = NaN;
    sigmoid_params = [];
    fit_quality = struct('success', false, 'rmse', NaN, 'rsquared', NaN);
    
    %% Step 1: Smooth the ACG using median filter
    filter_window_time = 0.83e-3; % 0.83 milliseconds in seconds
    filter_window_bins = round(filter_window_time / bin_size);
    if mod(filter_window_bins, 2) == 0
        filter_window_bins = filter_window_bins + 1; % Must be odd
    end
    
    acg_filtered = medfilt1(acg, filter_window_bins);
    
    %% Step 2: Find minimum (maximum value in 0-0.5ms window)
    % Note: "maximum value" is used to force the peak detection to be outside this window
    % This represents the end of the deepest part of the refractory period
    narrow_window = (bin_centers >= 0) & (bin_centers <= 0.5e-3);
    if ~any(narrow_window)
        warning('No bins in 0-0.5ms window');
        return;
    end
    
    min_value = max(acg_filtered(narrow_window));
    min_time = bin_centers(find(acg_filtered(narrow_window) == min_value, 1, 'last'));
    
    %% Step 3: Find peaks and select maximum
    [peak_values, peak_locs] = findpeaks(acg_filtered);
    peak_times = bin_centers(peak_locs);
    
    % Define threshold
    acg_range = max(acg_filtered) - min(acg_filtered);
    threshold = 0.1 * acg_range;
    
    % Keep peaks above threshold and above minimum
    valid_peaks = (peak_values > threshold) & (peak_values > min_value);
    
    if ~any(valid_peaks)
        warning('No valid peaks found');
        return;
    end
    
    valid_peak_times = peak_times(valid_peaks);
    valid_peak_values = peak_values(valid_peaks);
    
    % Select peak closest to 0 (earliest peak in recovery)
    [~, closest_idx] = min(abs(valid_peak_times));
    max_time = valid_peak_times(closest_idx);
    max_value = valid_peak_values(closest_idx);
    
    %% Step 4: Truncate ACG to fit window
    fit_window = (bin_centers >= min_time) & (bin_centers <= max_time);
    
    if sum(fit_window) < 3
        warning('Insufficient data points for fitting');
        return;
    end
    
    x_fit = bin_centers(fit_window);
    y_fit = acg_filtered(fit_window);
    
    %% Step 5: Fit sigmoid function
    % 4-parameter sigmoid: S(x) = ymin + (ymax - ymin) / (1 + exp(-k*(x - x0)))
    % This represents recovery from refractory period
    
    % Initial parameter estimates
    ymin_init = min_value;
    ymax_init = max_value;
    x0_init = mean([min_time, max_time]); % Midpoint of recovery
    k_init = 10 / (max_time - min_time); % Steepness estimate
    
    sigmoid_flexible = @(params, x) params(1) + (params(2) - params(1)) ./ ...
        (1 + exp(-params(3) * (x - params(4))));
    
    initial_params = [ymin_init, ymax_init, k_init, x0_init];
    
    % Bounds
    lb = [0, min_value, 0, min_time];
    ub = [max_value*2, max_value*2, Inf, max_time];
    
    % Fit options
    options = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxIterations', 1000);
    
    try
        [fitted_params, ~, residual, exitflag] = lsqcurvefit(sigmoid_flexible, ...
            initial_params, x_fit, y_fit, lb, ub, options);
        
        if exitflag <= 0
            warning('Sigmoid fit did not converge properly');
            return;
        end
        
        sigmoid_params = fitted_params;
        
        % Calculate goodness of fit
        y_pred = sigmoid_flexible(fitted_params, x_fit);
        fit_quality.rmse = sqrt(mean(residual.^2));
        SS_res = sum((y_fit - y_pred).^2);
        SS_tot = sum((y_fit - mean(y_fit)).^2);
        fit_quality.rsquared = 1 - SS_res/SS_tot;
        fit_quality.success = true;
        
    catch ME
        warning('Sigmoid fitting failed: %s', ME.message);
        return;
    end
    
    %% Step 6: Estimate refractory period
    % Find x such that S(x) = P_slope * (max(S) - min(S)) + min(S)
    % This is the time when neuron has recovered P_slope% of the way
    
    min_S = fitted_params(1);
    max_S = fitted_params(2);
    k = fitted_params(3);
    x0 = fitted_params(4);
    
    target_value = P_slope * (max_S - min_S) + min_S;
    
    % Solve for x analytically
    ratio = (max_S - min_S) / (target_value - min_S) - 1;
    
    if ratio <= 0
        warning('Invalid ratio in RP calculation');
        return;
    end
    
    RP_estimated_sec = x0 - log(ratio) / k;
    RP_estimated = RP_estimated_sec * 1000; % Convert to milliseconds
    
    %% Step 7: Check against first non-null bin
    first_nonzero_idx = find(acg > 0, 1, 'first');
    if ~isempty(first_nonzero_idx)
        first_nonzero_time = bin_centers(first_nonzero_idx) * 1000; % in ms
        if RP_estimated < first_nonzero_time
            RP_estimated = first_nonzero_time;
        end
    end
    
end