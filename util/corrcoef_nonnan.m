function [r,p] = corrcoef_nonnan(a,b)
if nargin == 1
    r = zeros(size(a,2)); p = zeros(size(a,2));
    for k = 1:size(a,2)
        for j = k:size(a,2)
            [nextr, nextp] = corrcoef_nonnan(a(:,k), a(:,j));
            r(k,j) = nextr(2); p(k,j) = nextp(2);
            r(j,k) = r(k,j); p(j,k) = p(k,j);
        end
    end
    return;
end
f = ~isnan(a) & ~isnan(b);
[r,p] = corrcoef(a(f),b(f));
