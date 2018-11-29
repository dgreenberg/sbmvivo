function o = ui_ison(h)
if ~all(ishandle(h))
    error('invalid handle');
end
o = false(size(h));
for j = 1:numel(h)
    switch get(h(j),'type')
        case 'uimenu'
            o(j) = strcmp(get(h(j),'checked'),'on');
        case 'uicontrol'
            if ismember(get(h(j),'style'),{'togglebutton','radiobutton','checkbox'})
                o(j) = logical(get(h(j),'value'));
            else
                error('invalid uicontrol type');
            end
        otherwise
            error('invalid object type');
    end
end