function [spikecounts_true_for_stats, spikecounts_est_for_stats, t_spikecounts_for_stats] = retrieve_spikecounts_for_stats(d, results)
%this function returns true and estimated spike counts with which to calculate statistics
%it tries to get them directly on the same time base from the data / results
%if this is not possible it will bin spike trains or resample existing counts'
[spikecounts_true_for_stats, spikecounts_est_for_stats, t_spikecounts_for_stats] = deal([]);

%get info about data
[~, it, dt, nstrue, ~, ~, sttrue, ~] = extract_data_fromoedatasegment(d);
twin_true = d.spiketimes_window;
truespiketimes_available = ~any(isnan(twin_true));

%get info about results
t_nsest = results.spikecount_times;
nsest = results.spikecounts;
stest = results.spiketimes;
twin_est = results.spiketimes_window;
dt_t_est = median(diff(t_nsest));

spikecounts_estimated = ~isempty(nsest);
spiketimes_estimated  = ~isempty(stest);

if spikecounts_estimated %the algorithm has output spike counts, so we use them
    
    if d.apcounts_present
        
        if all(ismember(t_nsest, it)) %estimated spike counts are available on the same time base as true spike counts
            
            t_spikecounts_for_stats = intersect(it, t_nsest);
            
            spikecounts_est_for_stats = nsest(ismember(t_nsest, it));
            
            spikecounts_true_for_stats = nstrue(ismember(it, t_nsest));
            
        elseif d.spiketimes_present %true spike train is available so we can bin it to get spike counts over the same time bins as estimated spike counts
            
            ii = find(t_nsest(:) >= twin_true(1) & t_nsest(:) <= twin_true(2));
            
            t_spikecounts_for_stats = t_nsest(ii);
            
            spikecounts_est_for_stats = nsest(ii);
            
            edges = [t_nsest(ii) - dt_t_est / 2; t_nsest(ii(end)) + dt_t_est / 2];
            counts = histc(sttrue, edges);
            spikecounts_true_for_stats = counts(1:end - 1);
            
        else %we have to re-bin the true and estimated spike counts to use the same time base
            
            %FIXME should use whatever is slower as the time base, right now we're using the true spikes (which also may have a longer time base, which is good)
            t_spikecounts_for_stats = it;
            
            spikecounts_true_for_stats = nstrue;
            ncounts = numel(spikecounts_true_for_stats);
            
            edges = [it(:) - dt / 2; it(end) + dt / 2];
            [~, bin_index] = histc(t_nsest, edges);
            spikecounts_est_for_stats = nan(ncounts, 1);
            
            for k = 1:ncounts
                
                spikecounts_est_for_stats(k) = mean(nsest(bin_index == k));
                
            end
            
        end
        
    end
    
elseif spiketimes_estimated %the algorithm has output spike trains but not spike counts
    
    if d.apcounts_present %use them
        
        ii = find(it(:) >= twin_est(1) & it(:) <= twin_est(2));
        edges = [it(ii) - dt / 2; it(ii(end)) + dt / 2];
        
        spikecounts_true_for_stats = nstrue(ii);
        
    elseif truespiketimes_available %bin true spike times around each scan time
        
        ii = find(it(:) >= twin_est(1) & it(:) <= twin_est(2) & it(:) >= twin_true(1) & it(:) <= twin_true(2));
        edges = [it(ii) - dt / 2; it(ii(end)) + dt / 2];
        
        counts = histc(sttrue, edges);
        spikecounts_true_for_stats = counts(1:end - 1);
        
    end
    
    t_spikecounts_for_stats = it(ii);
    
    %bin estimated spike train to give counts that match true spike counts in time
    counts = histc(stest, edges);
    spikecounts_est_for_stats = counts(1:end - 1);
    
end

%return column vectors:
[spikecounts_true_for_stats, spikecounts_est_for_stats, t_spikecounts_for_stats] = deal(spikecounts_true_for_stats(:), spikecounts_est_for_stats(:), t_spikecounts_for_stats(:));