function ui_setval(h,x, ndigits)
%setval(h,x)
%set string field(s) of GUI object(s)
%with handle(s) h to the string version of
%numbers x
if nargin < 3
    ndigits = 10;
end
if ~all(ishandle(h))
    error('invalid handle');
end
if numel(x) == 1
    x(1:numel(h)) = x;
end
for j=1:numel(h)
    ht = get(h(j),'type');
    if strcmp(ht,'uicontrol') || strcmp(ht,'text')
        if isempty(x) || isnan(x(j))
            set(h(j),'string','');
        else
            set(h(j),'string',num2str(x(j),ndigits));
        end
    else
        warning('handle is not to a text or uicontrol, ignoring');
    end
end