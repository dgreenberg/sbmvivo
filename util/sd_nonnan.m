function sd = sd_nonnan(x,flag,dim)
%sd = std_nonnan(x,dim)
%std over the given dimension 
%(or first nonsingleton one if none specified)
%skipping all NaN values
s = size(x);
if nargin < 2
    dim = min(find(s > 1));
    if isempty(dim)
        sd = nan + zeros(size(x));
        return;
    end
end
if nargin < 3
    flag = 0;
else
    assert(flag == 0 || flag == 1);
end

xnan = isnan(x);
x(find(xnan)) = 0; %set NaNs to zero

m = sum(x,dim) ./ sum(~xnan,dim);

% Avoid divide by zero for scalar case
if size(x,dim)==1, sd = zeros(size(x)); sd(isnan(x))=NaN; return, end

tile = ones(1,length(s));
tile(dim) = size(x,dim);

xc = x - repmat(m,tile);  % Remove mean
xc(find(xnan)) = 0; %set NaNs to zero
sd = sqrt(sum(xc.^2,dim) ./ (sum(~xnan,dim) - 1 + flag));
sd(find(~(sum(~xnan,dim) - 1 + flag))) = nan;