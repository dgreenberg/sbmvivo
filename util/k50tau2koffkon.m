function [koff, kon] = k50tau2koffkon(k50, tau)
%[koff, kon] = k50tau2koffkon(k50, tau)

n = numel(k50);

A = bsxfun(@power, k50(:), 1:n);
for i = 1:n
    
    A(i, 1:i - 1) = -A(i, 1:i - 1);
    
end

b = ones(n, 1);
beta = A \ b;
if any(beta < 0)
    
    beta = lsqnonneg(A, b);
    %FIXME: check values, issue a warning if a tolerance is exceeded
    
end
    
b = [1.0 beta'];
Ka = b(2:end) ./ b(1:end - 1);
kobs = 1.0 ./ tau;
koff = kobs ./ (1.0 + k50 .* Ka);
kon = Ka .* koff;