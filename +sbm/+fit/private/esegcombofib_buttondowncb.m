function esegcombofib_buttondowncb(figh, eventdata, fignames) %#ok<INUSL>
ax = findobj(figh, 'type', 'axes');
cp = get(ax, 'currentPoint');
cp = cp(1, 1:2);
XL = get(ax, 'xlim');
YL = get(ax, 'ylim');
if cp(1) < XL(1) || cp(1) > XL(2) || cp(2) < YL(1) || cp(2) > YL(2)
    
    return;
    
end
L = findobj(ax, 'type', 'line', 'marker', 'none');
x = get(L, 'xdata');
y = get(L, 'ydata');
[~, ii] = min((x - cp(1)) .^ 2 + (y - cp(2)) .^ 2);
delete(findobj(ax, 'color', 'g', 'marker', 'o'));
line(x(ii), y(ii), 'color','g','marker','o');
t = title(ax, fignames{ii});
set(t,'interpreter','none')