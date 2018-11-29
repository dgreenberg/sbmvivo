function x = tridiag_extrarowcol_solve(M, y)
% x = diag_extrarowcol_solve(M, y)
%
% solve M * x == y for x
%
% where M = [z v'; u D] and D is tridiagonal
%
% see https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
w = y(1);
d = y(2:end);
n = size(M, 1) - 1;
assert(ismatrix(M) && size(M, 2) == n + 1, 'M must be square');
D = M(2:end, 2:end);
a = [0; diag(D, -1)];
b = diag(D);
c = diag(D, 1);
assert(all(all(D == diag(a(2:end), -1) + diag(b) + diag(c, 1))), 'M(2:end, 2:end) must be tridiagonal');
z = M(1, 1);
u = M(2:end, 1);
v = M(1, 2:end);

wprime = w;
zprime = z;
cprime = nan(n - 1, 1);
[dprime, uprime] = deal(nan(n, 1));

cprime(1) = c(1) / b(1);
dprime(1) = d(1) / b(1);
uprime(1) = u(1) / b(1);

% we first make a forward pass to set a to zero and b to one
for i = 2:n
    
    % for loop iteration i, we add a multiple of row (i-1) to row i to set
    % a(i) to zero, then rescale row i to set b(i) to one. we also compute
    % the effect this has on c(i), d(i), and u(i)    
    
    denom = b(i) - a(i) * cprime(i - 1);
    
    dprime(i) = (d(i) - a(i) * dprime(i - 1)) / denom;
    uprime(i) = (u(i) - a(i) * uprime(i - 1)) / denom;
    
    if i < n
        
        cprime(i) = c(i) / denom;
        
    end
    
end

% we next make a backward pass to set c and v to zero
for i = n-1:-1:0
   
    if i > 0  % subtract cprime(i) * (row (i + 1)) of D from row i of D
        
        dprime(i) = dprime(i) - cprime(i) * dprime(i + 1);
        uprime(i) = uprime(i) - cprime(i) * uprime(i + 1);
        
    end
    
    % subtract v(i + 1) * (row (i + 1)) of D from row 1 of M
    zprime = zprime - v(i + 1) * uprime(i + 1);
    wprime = wprime - v(i + 1) * dprime(i + 1);
    
end

% rescale the first row of M to set zprime to one
wprime = wprime / zprime;

% finally we make a forward pass to set u to zero
dprime = dprime - wprime * uprime;

% [zprime vprime'; uprime diag(aprime(2:end), -1) + diag(bprime) + diag(cprime, 1)]
% is now the identity matrix

x = [wprime; dprime];