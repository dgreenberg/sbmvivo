function [y,baseline] = kinetic2dff_windowed_skipnan(r,p,w,blgsd)
if nargin < 4
    blgsd = 0;
end
y = nan + zeros(size(r));
if nargout > 1
    baseline = nan + zeros(size(r));
end
for k = 1:size(r,2)
    f = find(~isnan(r(:,k)));
    g = find(isnan(r(:,k)));
    s = unique([min(f); intersect(f,g + 1)]);
    e = unique([intersect(f,g - 1); max(f)]);
    for j = 1:length(s)
        if e(j) - s(j) < 1
            continue;
        elseif ~any(r(s(j):e(j),k))
            if nargout > 1
                baseline(s(j):e(j),k) = 0;
            end
            continue;
        end
        if nargout < 2
            y(s(j):e(j),k) = kinetic2dff_windowed(r(s(j):e(j),k),p,w,blgsd);
        else
            [y(s(j):e(j),k), baseline(s(j):e(j),k)] = kinetic2dff_windowed(r(s(j):e(j),k),p,w,blgsd);
        end
    end
end