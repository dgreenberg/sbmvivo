function [r,p] = spearman_rankcorr(x,y,alpha)
n = length(x);
rx = get_ranks(x);
ry = get_ranks(y);
d2 = (rx - ry).^2;
Tx = calc_T(x);
Ty = calc_T(y);
n3nd6 = (n^3 - n) / 6;
r = (n3nd6 - sum(d2) - Tx - Ty) / sqrt((n3nd6 - 2 * Tx) * (n3nd6 - 2 * Ty));
if nargout > 1
    if nargin < 3
        alpha = 0.01;
    end
    rvals = zeros(ceil(10 / alpha),1);
    for k = 1:size(rvals,1)
        rvals(k) = spearman(x(randperm(n)),y);
    end
    p = 1 - mean(abs(rvals) < abs(r));
end

function T = calc_T(a)
T = 0;
a = sort(a);
u = unique(a(find(diff(a) == 0)));
for j = 1:length(u)
    s = sum(a == u(j));
    T = T + s^3 - s;
end
T = T / 12;

function rx = get_ranks(x)
s = sortrows([x (1:length(x))'],1);
rx = zeros(length(x),1);
rx(s(:,2)) = (1:length(x))';
u = unique(s(find(diff(s(:,1)) == 0),1));
for j = 1:length(u)
    f = find(x == u(j));
    rx(f) = mean(rx(f));
end
