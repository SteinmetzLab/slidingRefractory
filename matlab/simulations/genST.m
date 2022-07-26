

function st = genST(rate, duration)
% function st = genST(rate, duration)
%
% generates a spike train with specified rate and duration. Rate and
% duration should have the same units. 

%%
mu = 1/rate;
n = rate*duration; 
isi = exprnd(mu, ceil(n*2), 1); % generate twice as many spikes as likely required

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

