function [k50, tau] = koffkon2k50tau(koff, kon)
%[k50, tau] = koffkon2k50tau_bs(koff, kon)
%convert from rate constants to K50s and time constants
n = numel(kon);
Ka = kon ./ koff;
beta = [1 cumprod(Ka)];
% given a set of binding constants beta, determine the free ligand
% concentrations gamma for which each binding step is half completed.
% works by root finding.
k50 = nan(size(koff));
for i = 1:n
    
    s = (-1) .^ ((0:n) < i);
    Q = beta .* s;
    g = roots(Q(end:-1:1));
    ok = g > 0 & imag(g) == 0;
    assert(sum(ok) == 1, 'failed to find unique positive real root');
    k50(i) = g(ok);
    
end
kobs = koff .* (1.0 + k50 .* Ka);
tau = 1.0 ./ kobs;
    
    
    