function x = listboxchoice(h)
%function listchoice(h)
%gets the chosen string
%from listbox or popup menu
%with handle h
c=get(h,'string');
if iscell(c)
    x = c{get(h,'value')};
else
    x = c;
end
x(double(x) == 13) = []; %remove carraige returns