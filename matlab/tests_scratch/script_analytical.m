

%% test a closed form expression for probability of passing
clc; clear all;

PTrue = 0.05; % proportion true contamination
PLimit = 0.1; % max acceptable proportion contamination
confLimit = 0.9; % confidence that we find acceptable
R = 0.9; % firing rate
% r = 0.002; % refractory period
D = 3*3600; % recording duration

rpBinSize = 1/30000;
rValues = 0.5/1000:rpBinSize:10/1000; % refractory periods to check over

verbose = true;

testR = logspace(-3,1.5,100); 

allPropPass = zeros(numel(testR), numel(rValues)); 

for Ridx = 1:numel(testR)
    R = testR(Ridx); 

    S = R*D; % spike count
    
    for ridx = 1:numel(rValues)
        r = rValues(ridx); 
    
        violLimit = PLimit*R*r*2*S; % mean count of a poisson process at the limit of acceptable contamination
        violTrue = PTrue*R*r*2*S; % mean count of poisson process at the true contamination

        % what is the value of the limit contamination's CDF that would exceed the
        % acceptable confidence? I.e. when does the limit CDF = 10%? this is the
        % max acceptable. 
        k = poissinv(1-confLimit, violLimit);     

        % now what is the value of the true CDF just below that value? that's the
        % proportion of simulations that should pass
        allPropPass(Ridx, ridx) = poisscdf(k-1, violTrue); % k can be negative, in which case it returns zero
    end
    
    % now, if the test at every refractory period was completely
    % independent, this would be the probability of passing
    propPass(Ridx) = 1-prod(1-allPropPass(Ridx,:)); 
    
    if verbose
        fprintf(1, 'rate = %.2f, spk count = %.0f, limit viol = %.2f, true mean viol = %.2f, max count = %d, pct pass=%.2f\n', ...
            R, S, violLimit, violTrue, k, propPass(Ridx)*100);
    end
end
figure; semilogx(testR, propPass*100, '.-');
xlabel('firing rate (sp/s)'); 
ylabel('percent passing'); 

% there's a funny effect where this function is non-monotonic due to the
% discrete nature of the spike counts: if you're allowed 1 contaminating
% spike (for instance) then any increases in spike rate up to the point
% where you're allowed 2 contaminating spikes actually result in decreased
% probability to pass, since you have more contamination (in raw numbers),
% so less probability of getting count of zero or 1. I think that's a
% correct outcome, as weird as it seems. shouldn't make a difference in any
% practical sense. 

% now let's put the above into a function and test it on a range of values.

%% tests

% iterate over:
% pct contamination (rows)
% true RP (columns)
% 



