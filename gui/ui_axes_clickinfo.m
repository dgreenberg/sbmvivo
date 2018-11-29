function [clickax, clickdatapos, Lclick, Rclick, Oclick] = ui_axes_clickinfo(figh)
Lclick = strcmpi(get(figh,'selectiontype'),'normal');
Rclick = strcmpi(get(figh,'selectiontype'),'alt');
Oclick = strcmpi(get(figh,'selectiontype'),'open');
clickax = findobj(figh,'type','axes');
clickax = clickax(find(ui_hasmouse(clickax), 1));
if isempty(clickax)
    clickdatapos = [NaN NaN];
else
    cp = get(clickax,'currentpoint');
    clickdatapos = cp(1,1:2);
end