function oedb_compareresults(oedb, fileinfo, dataset, neuron, segment, XL)
resultspresent = find(~cellfun(@isempty, {oedb.stats.bydataset(:, dataset).corr_spikecounts}));
if isempty(resultspresent), return; end
[f,it,dt,ns,~,~,st,~,~,stwin] = oedb_fetch_timeseries(oedb, fileinfo, dataset, neuron, segment);
if ~exist('XL', 'var') || isempty(XL)
    
    XL = it([1 end]) + [-1 1] * dt / 2;
    
end

nax = numel(resultspresent) + 2;
fh = figure('renderer', 'painters');
ax(1) = subplot(nax, 1, 1, 'parent', fh);
plot(it, f, 'k');
ax(2) = subplot(nax, 1, 2, 'parent', fh);
if ~any(isnan(stwin))
    
    plotspiketrain(st, ax(2));
    set(ax(2), 'ytick', []);

else
    
    plot(it, ns, 'parent', ax(2));
    
end
ylabel('true APs');

for a = 1:numel(resultspresent)
   
    algind = resultspresent(a);
    results = fetch_results(oedb, fileinfo, dataset, neuron, segment, algind);
    
    ax(a + 2) = subplot(nax, 1, a + 2, 'parent', fh);
    
    if ~any(isnan(results.spiketimes_window))
        
        plotspiketrain(results.spiketimes, ax(a + 2));
        set(ax(a + 2), 'ytick', []);
        
    else
        
        line(results.spikecount_times, results.spikecounts, 'color', 'k', 'parent', ax(a + 2));
        
    end
    
    ylabel(oedb.algnames{algind});
    
    if a < numel(resultspresent)
       
        set(ax(a + 2), 'xtick', []);
        
    end
    
end
linkaxes(ax, 'x');
set(ax, 'xlim', XL, 'box', 'off', 'tickdir', 'out');