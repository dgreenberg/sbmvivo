function plot_f_and_spikes(it, f, st, ns, ax, f_removed, showremoved)
if ~exist('ax', 'var') || isempty(ax)
   
    figure;
    ax(1) = subplot(2, 1, 1);
    ax(2) = subplot(2, 1, 2);
    
end
if exist('f_removed', 'var') && ~isempty(f_removed)
    if ~exist('showremoved', 'var') || showremoved
        
        if size(f_removed, 2) > 1
            
            line(it, bsxfun(@plus,f,f_removed), 'parent', ax(1));
            
        else
            
            line(it, f + f_removed, 'parent', ax(1), 'color', [1 0.5 0]);
            
        end
        
    end
end
line(it, f, 'parent', ax(1), 'color', [0 .7 0],'marker','.');
if ~isempty(st) %plot true (electrically detected) spike train
    
    plotspiketrain(st, ax(2), 2.1, 0.9);
    
elseif ~isempty(ns) %if true spike train is unavailable, plot true spike counts instead
    
    line(it, ns, 'color', 'k', 'parent', ax(2));
    
end