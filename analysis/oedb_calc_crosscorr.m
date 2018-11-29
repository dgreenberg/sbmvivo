function [crosscorr, shift_s] = oedb_calc_crosscorr(maxshift_s, shiftstep_s, sigma_s, oedb, fileinfo, algind, dataset, neuron, bin_size_s)
%crosscorr = oedb_calc_crosscorr(maxshift_s, shiftstep_s, sigma_s, oedb, fileinfo, algind, dataset, neuron, bin_size_s)
if ~exist('bin_size_s', 'var') || isempty(bin_size_s)
    bin_size_s = 1e-3;
end
oerec = fetch_neurondata(oedb, fileinfo, dataset, neuron);
nseg = numel(oerec.data);
[spiketimes_true_for_stats, spiketimes_est_for_stats, weights_true_for_stats, weights_est_for_stats] = deal(cell(1, nseg));
twin_for_stats = nan(nseg, 2);

%retrieve true and predicted spike times
for si = 1:nseg
    
    data = oerec.data(si);
    results = fetch_results(oedb, fileinfo, dataset, neuron, si, algind);
    
    %retrieve true and estimated spike times / window:
    [spiketimes_true_for_stats{si}, spiketimes_est_for_stats{si}, weights_true_for_stats{si}, weights_est_for_stats{si}, twin_for_stats(si, :)] = retrieve_spiketimes_for_stats(data, results);
    %note that only spiketimes within the relevant windows are returned
    
end

maxsteps = ceil(maxshift_s / shiftstep_s);
shift_s = (-maxsteps:maxsteps)' * shiftstep_s;
crosscorr = nan(2 * maxsteps + 1, 1);
for j = 1:numel(crosscorr)
    
    sa = spiketimes_true_for_stats;
    sb = cellfun(@(v) v + shift_s(j), spiketimes_est_for_stats, 'uniformoutput', false);
    twin_shifted = [max(twin_for_stats(:, 1), twin_for_stats(:, 1) + shift_s(j)), min(twin_for_stats(:, 2), twin_for_stats(:, 2) + shift_s(j))]; % contract the window
    [wa, wb, ns_a, ns_b, timebase] = deal(cell(1, nseg));
    for k = 1:nseg
        
        ok_a = sa{k} >= twin_shifted(k, 1) & sa{k} <= twin_shifted(k, 2);
        sa{k} = sa{k}(ok_a);
        wa{k} = weights_true_for_stats{k}(ok_a);
        ok_b = sb{k} >= twin_shifted(k, 1) & sb{k} <= twin_shifted(k, 2);
        sb{k} = sb{k}(ok_b);
        wb{k} = weights_est_for_stats{k}(ok_b);
        
        bin_edges = [twin_shifted(k, 1):bin_size_s:twin_shifted(k, 2) + bin_size_s];
        timebase{k} = bin_edges(1:end-1) + bin_size_s / 2;
        
        [ns_a{k}, ns_b{k}] = deal(zeros(numel(bin_edges) - 1, 1));
        [~, bin_index_a] = histc(sa{k}, bin_edges);
        [~, bin_index_b] = histc(sb{k}, bin_edges);
        
        for u = 1:numel(sa{k})
            
            ns_a{k}(bin_index_a(u)) = ns_a{k}(bin_index_a(u)) + wa{k}(u);
            
        end
        for u = 1:numel(sb{k})
            
            ns_b{k}(bin_index_b(u)) = ns_b{k}(bin_index_b(u)) + wb{k}(u);
            
        end
    end
    
    [~, ~, ~, ~, crosscorr(j)] = oedb_masked_filter(ns_a, ns_b, sigma_s, timebase);
    
end