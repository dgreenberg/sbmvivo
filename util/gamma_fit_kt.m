function [kt, LL] = gamma_fit_kt(x, options)
% [kt, LL] = gamma_fit_kt(x)
% see https://en.wikipedia.org/wiki/Gamma_distribution#Maximum_likelihood_estimation
assert(isnumeric(x), 'x must be numeric');
x = double(x(:));
assert(numel(x) > 0, 'no data');
assert(all(x > 0) && all(imag(x) == 0) && all(isfinite(x)), 'values must be positive');
if ~exist('options','var') || isempty(options)
   
    options = optimset('Display','none','maxiter',1000,'tolx',1e-8);
    
end

N = numel(x);
sum_log_x = sum(log(x));
log_mean_x = log(mean(x));

s = log_mean_x - sum_log_x / N;
k0 = (3 - s + sqrt((s - 3) ^ 2 + 24 * s)) / (12 * s);
if k0 < 0
    
    k0 = 1;
    
end

ofunc = @(k) -LL_autot(k, sum_log_x, log_mean_x, N);
[k, negLL] = fminbnd(ofunc, 0.1 * k0, 2 * k0, options);

t = mean(x) / k;
LL = -negLL;
kt = [k t];


function LL = LL_autot(k, sum_log_x, log_mean_x, N)
LL = (k - 1) * sum_log_x ...
    - N * k ...
    - N * k * (log_mean_x - log(k)) ...
    - N * gammaln(k);