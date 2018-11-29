function ui_toggleon(h)
if ~all(ishandle(h))
    error('invalid handle');
end
wason = ui_ison(h);
for j = 1:numel(h)
    switch get(h(j),'type')
        case 'uimenu'
            if wason(j)
                set(h(j), 'checked', 'off');
            else
                set(h(j), 'checked', 'on');
            end            
        case 'uicontrol'
            if ismember(get(h(j),'style'),{'togglebutton','radiobutton','checkbox'})
                set(h(j), 'value', double(~wason(j)));                
            else
                error('invalid uicontrol type');
            end
        otherwise
            error('invalid object type');
    end
end