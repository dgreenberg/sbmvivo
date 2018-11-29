function m = ui_hasmouse(a)
%function m = hasmouse(a)
%returns true if the currentpoint of 2D axes with handle a falls inside its
%axis limits (inclusive), false otherwise or if axis not visible and has no child objects
if numel(a) > 1
    m = nan(size(a));
    for u = 1:numel(a)
        m(u) = ui_hasmouse(a(u));
    end
    return;
end
if strcmp(get(a,'visible'),'off') && isempty(get(a,'children'))
    m = false;
    return;
end
f = get(a,'parent'); %figure handle
ou = get(f,'units');
set(f,'units',get(0,'units'));
p = get(f,'position');
c0 = get(0,'PointerLocation');
set(f,'units',ou);
% if c0(1) < p(1) || c0(1) > p(1) + p(3) || c0(2) < p(2) || c0(2) > p(2) + p(4)
%     m = false;
%     return;
% end
c = get(a,'currentpoint');
x = c(1,1);
y = c(1,2);
XL = get(a,'xlim');
YL=get(a,'ylim');
if (x < XL(1)) || (x > XL(2)) || (y < YL(1)) || (y > YL(2))
    m = false;
else
    m = true;
end