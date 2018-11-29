function m = mean_nonnan(x,dim)
%m = mean_nonnan(x,dim)
%mean over the given dimension 
%(or first nonsingleton one if none specified)
%skipping all NaN values
if nargin < 2
    
    s = size(x);
    dim = find(s > 1, 1 );
    if isempty(dim)
        if isempty(x)
            m = [];
        else
            m = x;
        end
        return;
    end
    
end

xnan = isnan(x);
x(xnan) = 0; %set NaNs to zero

m = sum(x, dim) ./ sum(~xnan, dim);