function n = partition_vector(v,r)
v = reshape(v,[],1);
n = zeros(size(v));
fv = find(v == 1);
if isempty(fv)
    return;
end
dfv = diff(fv);
g = find(dfv > r);
s = [fv(1); fv(g + 1)];
e = [fv(g); fv(end)];
for k = 1:length(s)
    n(s(k):min(end,e(k) + r)) = k;
end