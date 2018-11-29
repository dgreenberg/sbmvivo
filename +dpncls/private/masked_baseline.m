function b = masked_baseline(y, mask, ws, s)
kernel = exp( -[(-ws:ws).^2'/(2*s^2)] ) / (s * sqrt(2*pi));

n = length(y);
y = [zeros(ws,1); y; zeros(ws,1)];
mask = [zeros(ws,1); mask; zeros(ws,1)];

indmat = repmat(ws + 1:ws + n, 2 * ws + 1, 1) + repmat((-ws:ws)',1,n);

im = mask(indmat);
iy = y(indmat);
ik = repmat(kernel,1,n);

iy(~im) = 0;
ik(~im) = 0;

warning off MATLAB:divideByZero
b = (sum(iy .* ik) ./ sum(ik))';

mask = mask(1 + ws:end-ws);
y = y(1 + ws:end - ws);

fm = find(mask);
for k = find(isnan(b) & ~isnan(y))'
    f = find(fm > k);
    g = find(fm < k);
    z = [];
    if ~isempty(f)
        z = f(1:min(end,3));
    end
    if ~isempty(g)
        z = [z; g(max(1,end-2):end)];
    end
    if isempty(z)
        b(k) = nan;
    else
        b(k) = mean_nonnan(y(fm(z)));
    end
end
warning on MATLAB:divideByZero

