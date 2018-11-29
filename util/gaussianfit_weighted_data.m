function [mu, sigma] = gaussianfit_weighted_data(x,w)
%[mu, sigma] = gaussianfit_weighted_data(x)
sw = sum(w);
mu = w' * x / sw;
r = x - mu;
resid_ssq = w' * (r .^ 2);
sigma = sqrt(resid_ssq / sw);