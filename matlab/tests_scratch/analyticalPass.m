
function propPass = analyticalPass(...
    PTrue,PLimit,confLimit,R,r,D,verbose)
% input arguments are defined as follows, with example values
% PTrue = 0.05; % proportion true contamination
% PLimit = 0.1; % max acceptable proportion contamination
% confLimit = 0.9; % confidence that we find acceptable
% R = 0.9; % firing rate
% r = 0.002; % refractory period
% D = 3*3600; % recording duration
% verbose = true;
if nargin<7
    verbose = false;
end

S = R*D; % spike count
violLimit = PLimit*R*r*2*S; % mean count of a poisson process at the limit of acceptable contamination
violTrue = PTrue*R*r*2*S; % mean count of poisson process at the true contamination

% what is the value of the limit contamination's CDF that would exceed the
% acceptable confidence? I.e. when does the limit CDF = 10%? this is the
% max acceptable.
k = poissinv(1-confLimit, violLimit);

% now what is the value of the true CDF just below that value? that's the
% proportion of simulations that should pass
propPass = poisscdf(k-1, violTrue); % k can be negative, in which case it returns zero

if verbose
    fprintf(1, 'rate = %.2f, spk count = %.0f, limit viol = %.2f, true mean viol = %.2f, max count = %d, pct pass=%.2f\n', ...
        R, S, violLimit, violTrue, k, propPass(Ridx)*100);
end