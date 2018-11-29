function [w, log_w_norm, neff, log_sum_raw_weights] = normalize_weights(log_w)
%[w, log_w, neff, log_sum_raw_weights] = normalize_weights(log_w)
%given logarithms of raw weights log_w, this function uses careful
%calculations to computes the following quantities while seeking to
%minimize numerical instability by working with sums of very large or small
%numbers in the linear domain:
%
%w                   - normalized weights in the linear domain
%log_w_norm          - logarithms of normalized weights
%neff                - effective sample size = 1 / sum(w .^ 2)
%log_sum_raw_weights - sum of log_w, the standard sequential monte carlo estimate of P(Y(t) | Y(1:t - 1))
log_w_max           = max(log_w);
w                   = exp(log_w - log_w_max);
sw                  = sum(w);
lsw                 = log(sw);
log_w_norm          = log_w - log_w_max - lsw;
w                   = w / sw;
neff                = 1 / sum(w .^ 2);
log_sum_raw_weights = lsw + log_w_max;