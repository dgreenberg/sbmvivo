function spikecounts_split = getspikecounts(stsplit, spiketshift, bin_edges)
for u = 1:numel(stsplit)
    if ~(all(stsplit{u} + spiketshift >= bin_edges{u}(1)) && all(stsplit{u} + spiketshift <= bin_edges{u}(end)))
        %warning('some spikes lost due to clipping');
    end
    spikecounts_split{u} = histc(stsplit{u} + spiketshift, bin_edges{u});
    spikecounts_split{u}(end) = []; %last bin is only the counts matching the final edge precisely
end