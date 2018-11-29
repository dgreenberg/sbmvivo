function [hline, hpatch] = msdplot(x, y, yplus, yminus, parent, linecolor, patchcolor)

x = x(:);
y = y(:);
yplus = yplus(:);
yminus = yminus(:);

if exist('parent', 'var')
    
    newfig = false;
    
else
    
    newfig = true;
    fig = figure;
    parent = axes('parent', fig);
    
end

if ~exist('linecolor', 'var')
    
    linecolor = [0 0 0];
    
end
if ~exist('patchcolor', 'var')
    
    patchcolor = [.5 .5 .5];
    
end

ok = isfinite(x) & isfinite(yplus) & isfinite(yminus);
blockstarts = find(diff([false; ok]) == 1);
blockends = find(diff([ok; false]) == -1);
hpatch = nan(1, numel(blockstarts));
for j = 1:numel(blockstarts)
    
    next_x = x(blockstarts(j):blockends(j));
    next_yp = yplus(blockstarts(j):blockends(j));
    next_ym = yminus(blockstarts(j):blockends(j));
    hpatch(j) = patch( ...
        [next_x; flipud(next_x); next_x(1)], ...
        [next_ym; flipud(next_yp); next_ym(1)], ...
        patchcolor, 'parent', parent, 'edgecolor', 'none');
    
end

hline = line(x, y, 'parent', parent, 'color', linecolor, 'linewidth', 2);

if newfig
    
    set(parent, 'xlimmode', 'auto', 'ylimmode', 'auto');
    
end