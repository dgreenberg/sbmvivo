function y = rms(x, dim)
if ~exist('dim', 'var')
    
    dim = find(size(x) > 1, 1);
    if isempty(dim)
        
        dim = 1;
        
    end
    
end
y = sqrt(mean(abs(x).^2, dim));