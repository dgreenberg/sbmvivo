function L = plotspiketrain(spiketimes, axh, y, h)
%L = plotspiketrain(spiketimes, axh, y, h)
if ~exist('axh','var') || isempty(axh)
    
    f = figure;
    axh = axes('parent', f);
    
end

if ~exist('y', 'var') || isempty(y)
    
    y = 0;
    
end

if ~exist('h', 'var') || isempty(h)
   
    h = 1;
    
end

xd = reshape([1; 1; nan] * spiketimes(:)', [], 1);
yd = reshape([y; y + h; nan] * ones(1, numel(spiketimes)), [], 1);
L = line(xd, yd, 'color', 'k', 'parent', axh);