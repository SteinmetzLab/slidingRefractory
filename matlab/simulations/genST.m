

function st = genST(rate, duration, varargin)
% function st = genST(rate, duration)
%
% generates a spike train with specified rate and duration. Rate and
% duration should have the same units. 

if nargin<3
    refPeriod = 0;
else
    refPeriod = varargin{1};
end

%%
% mu = 1/rate;
% n = rate*duration; 
% isi = exprnd(mu, ceil(n*2), 1); % generate twice as many spikes as likely required


rSim = rate/(1-refPeriod*rate);

mu = 1/rSim;

n = rate*duration; 

isi = refPeriod+exprnd(mu, ceil(n*2), 1); 

while sum(isi)<duration
    % this part will be extremely slow if it needs to be invoked, but in
    % general it should not be 
    isi(end+1) = exprnd(mu);
end

st = cumsum(isi); 
st = st(1:find(st<duration,1,'last'));

return;

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
