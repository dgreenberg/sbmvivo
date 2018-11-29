function [spiketimes_true_for_stats, spiketimes_est_for_stats, weights_true_for_stats, weights_est_for_stats, twin_for_stats] = retrieve_spiketimes_for_stats(d, results)
%this function returns true and estimated spike counts with which to calculate statistics
%it's a bit simpler than calculating spike counts since we don't need a common timebase, just a common time window
[spiketimes_true_for_stats, spiketimes_est_for_stats, weights_true_for_stats, weights_est_for_stats] = deal([]);
twin_for_stats = nan(1, 2);

%get info about data:
[~, it, ~, nstrue, ~, ~, sttrue, ~] = extract_data_fromoedatasegment(d);
twin_true = d.spiketimes_window;

%get info about results:
t_nsest = results.spikecount_times;
nsest = results.spikecounts;
stest = results.spiketimes;
twin_est = results.spiketimes_window;

spikecounts_estimated = ~isempty(nsest);
spiketimes_estimated  = ~any(isnan(twin_est));

if ~(spikecounts_estimated || spiketimes_estimated) || ~(d.spiketimes_present || d.apcounts_present)
    
    return; %need both data and results to do stats
    
end

if spiketimes_estimated
    
    spiketimes_est_for_stats = stest;
    weights_est_for_stats = ones(size(stest));
    twin_for_stats(1) = max(twin_for_stats(1), twin_est(1));
    twin_for_stats(2) = min(twin_for_stats(2), twin_est(2));
    
elseif spikecounts_estimated
    
    [spiketimes_est_for_stats, weights_est_for_stats, twin_sc] = convert_spikecounts_to_weighted_spiketrain(nsest, t_nsest);
    twin_for_stats(1) = max(twin_for_stats(1), twin_sc(1));
    twin_for_stats(2) = min(twin_for_stats(2), twin_sc(2));
    
end

if d.spiketimes_present
    
    spiketimes_true_for_stats = sttrue;
    weights_true_for_stats = ones(size(sttrue));
    twin_for_stats(1) = max(twin_for_stats(1), twin_true(1));
    twin_for_stats(2) = min(twin_for_stats(2), twin_true(2));
    
elseif d.apcounts_present
    
    [spiketimes_true_for_stats, weights_true_for_stats, twin_sc] = convert_spikecounts_to_weighted_spiketrain(nstrue, it);
    twin_for_stats(1) = max(twin_for_stats(1), twin_sc(1));
    twin_for_stats(2) = min(twin_for_stats(2), twin_sc(2));
   
end

%filter out estimated / true spike times outside the relevant time windows
ii = spiketimes_est_for_stats >= twin_for_stats(1) & spiketimes_est_for_stats <= twin_for_stats(2);
spiketimes_est_for_stats = spiketimes_est_for_stats(ii);
weights_est_for_stats = weights_est_for_stats(ii);

ii = spiketimes_true_for_stats >= twin_for_stats(1) & spiketimes_true_for_stats <= twin_for_stats(2);
spiketimes_true_for_stats = spiketimes_true_for_stats(ii);
weights_true_for_stats = weights_true_for_stats(ii);