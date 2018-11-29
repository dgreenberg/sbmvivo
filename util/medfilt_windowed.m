function y = medfilt_windowed(x, n)
%y = medfilt_windowed(x, n)
%n is full window size, and will be rounded up if not odd
%padding is done using the boundary elements
n = max(n, 3);
if mod(n, 2) == 0
   
    n = n + 1;
    
end
ws = (n - 1) / 2;
ii = bsxfun(@plus, (-ws:ws)', 1:numel(x));
ii = max(1, min(numel(x), ii));
xx = sort(x(ii));
y = reshape(xx(ws + 1, :), size(x));