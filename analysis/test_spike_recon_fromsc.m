function [npairs, nspikes_true, nspikes_est, Tseg] = test_spike_recon_fromsc(spikecounts_est_for_stats, spikecounts_true_for_stats, t_spikecounts_for_stats, tol)
nseg = numel(spikecounts_est_for_stats);

%first, we construct weighted spike trains from spike counts
[st_fromsc_est, w_fromsc_est, st_fromsc_true, w_fromsc_true] = deal(cell(1, nseg));
Tseg = nan(1, nseg);

for s = 1:nseg
    
    if isempty(spikecounts_est_for_stats{s}) || isempty(spikecounts_true_for_stats{s})
        
        continue;
        
    end
    
    [st_fromsc_est{s}, w_fromsc_est{s}] = convert_spikecounts_to_weighted_spiketrain(spikecounts_est_for_stats{s}, t_spikecounts_for_stats{s});
    [st_fromsc_true{s}, w_fromsc_true{s}] = convert_spikecounts_to_weighted_spiketrain(spikecounts_true_for_stats{s}, t_spikecounts_for_stats{s});
    
    dt = median(diff(t_spikecounts_for_stats{s}));
    Tseg(s) = max(t_spikecounts_for_stats{s}) - min(t_spikecounts_for_stats{s}) + dt;
    
end
%next, we calculate dr/fp as for spike trains
[npairs, nspikes_true, nspikes_est] = test_spike_recon(st_fromsc_true, st_fromsc_est, tol, w_fromsc_true, w_fromsc_est);