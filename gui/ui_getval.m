function x = ui_getval(h)
%x = getval(h)
%get numeric value(s) from string field
%of GUI object(s) with handle(s) h
if ~all(ishandle(h))
    error('invalid handle');
end
x = nan + zeros(size(h));
for j=1:numel(h)
    ht = get(h(j),'type');
    if strcmp(ht,'uicontrol') || strcmp(ht,'text')
        v = str2num(get(h(j),'string'));
        if ~isempty(v)
            x(j) = v;        
        end
   else
        warning('handle is not to a text or uicontrol, ignoring');
   end
end