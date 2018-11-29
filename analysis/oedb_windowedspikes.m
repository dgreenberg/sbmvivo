function [mean_estspikes, sd_estspikes] = oedb_windowedspikes(oedb, fileinfo, dsind, algind, winsize, edgebuffer, overlap_frac)
if ~exist('winsize', 'var'), winsize = 1; end
if ~exist('edgebuffer', 'var'), edgebuffer = [1; 1]; end % amount of trace to exlude at start/stop of recording
if ~exist('overlap_frac', 'var'), overlap_frac = 0; end %how much each window overlaps with the next
assert(overlap_frac >= 0 && overlap_frac < 1, 'overlap fraction must be nonnegative and less than 1');

n_overlap_steps = max(1, round(1 / (1 - overlap_frac)));

[mean_estspikes, sd_estspikes] = deal([]); %as a function of true spikes. first row is for zero true spikes in the window

for neuronind = 1:oedb.nneurons(dsind)
    
    oerec = fetch_neurondata(oedb, fileinfo, dsind, neuronind);
    
    true_and_est_counts = zeros(0, 2);
    
    for segmentind = 1:oedb.nsegments{dsind}(neuronind)
        
        data = oerec.data(segmentind);
        results = fetch_results(oedb, fileinfo, dsind, neuronind, segmentind, algind);
        
        [spiketimes_true_for_stats, spiketimes_est_for_stats, weights_true_for_stats, weights_est_for_stats, twin_for_stats] = ...
            retrieve_spiketimes_for_stats(data, results);
        
        edges = twin_for_stats(1) + edgebuffer(1):winsize:twin_for_stats(2) - edgebuffer(2);
        
        for v = 1:n_overlap_steps
            
            next_true_and_est_counts = nan(numel(edges) - 1, 2);
            
            edges_offset = winsize * (v - 1) / n_overlap_steps;
            next_edges = edges + edges_offset;
            next_edges(next_edges > twin_for_stats(2) - edgebuffer(2)) = [];
            
            [~, ii_true] = histc(spiketimes_true_for_stats, next_edges);
            [~, ii_est] = histc(spiketimes_est_for_stats, next_edges);
            
            for u = 1:numel(next_edges) - 1
                
                next_true_and_est_counts(u, 1) = sum(weights_true_for_stats(ii_true == u));
                next_true_and_est_counts(u, 2) = sum(weights_est_for_stats(ii_est == u));
                
            end
            
            true_and_est_counts = [true_and_est_counts; next_true_and_est_counts];
            
        end
        
    end
    
    max_truespikes = max(true_and_est_counts(:, 1));
    %grow the arrays:
    if neuronind > 1
        mean_estspikes(end + 1:max_truespikes + 1, :) = nan; %assign nan values if no previous neurons had this many true spikes
        sd_estspikes(end + 1:max_truespikes + 1, :) = nan; %assign nan values if no previous neurons had this many true spikes
    end
    mean_estspikes(:, end + 1) = nan; %add a new neuron
    
    sd_estspikes(:, end + 1) = nan; %add a new neuron
    
    for k = 0:max_truespikes
        
        y = true_and_est_counts(true_and_est_counts(:, 1) == k, 2);
        assert(~any(isnan(y)), 'unexpected NaN value')
        mean_estspikes(k + 1, end) = mean(y);
        sd_estspikes(k + 1, end) = std(y, 1);
        
    end
    
end