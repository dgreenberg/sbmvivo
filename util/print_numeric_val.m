function print_numeric_val(v, ndigits, nchar, nchar_multi)
if ~exist('nchar_multi','var')
    nchar_multi = 4;
end
if size(v,2) > 1
    for q = 1:size(v,2)
         print_numeric_val(v(:,q), ndigits, nchar_multi);
    end
    return;
end
if ischar(v)
    s = v;
elseif size(v,1) == 4
    s = [num2str(v(1),3) ' ' num2str(v(2),3) ' ' num2str(v(3),3) ' ' num2str(v(4),3)];
elseif size(v,1) > 1
    s = 'vector';
else
    s = num2str(v, ndigits);
end
s = [repmat(' ', 1, max(0, nchar - numel(s))) s];
fprintf(' %s',s); %note the extra space