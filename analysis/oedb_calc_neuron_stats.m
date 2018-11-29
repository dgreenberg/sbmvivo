function oedb = oedb_calc_neuron_stats(oedb, fileinfo, algind, dataset, neuron, segments_to_use)
%we have for each bit of data either
%A)various spike times / windows, or
%B)spikecounts / times
%we have this for both true and estimated spike counts
%for each segment of this neuron, we're going to try to retrieve both.
%when one is absent, we convert the other

%initialize
if ~exist('segments_to_use', 'var') || isempty(segments_to_use), segments_to_use = 1:oedb.nsegments{dataset}(neuron); end
nseg = numel(segments_to_use);

[spikecounts_true_for_stats, spikecounts_est_for_stats, t_spikecounts_for_stats, ...
    spiketimes_true_for_stats, spiketimes_est_for_stats, weights_true_for_stats, weights_est_for_stats] = deal(cell(1, nseg));
twin_for_stats = nan(nseg, 2);

oerec = fetch_neurondata(oedb, fileinfo, dataset, neuron);

%retrieve true and predicted spike counts/times
for si = 1:nseg
    
    s = segments_to_use(si);
    data = oerec.data(s);
    results = fetch_results(oedb, fileinfo, dataset, neuron, s, algind);
    
    %retrieve true and estimated spike counts / time base:
    [spikecounts_true_for_stats{si}, spikecounts_est_for_stats{si}, t_spikecounts_for_stats{si}] = retrieve_spikecounts_for_stats(data, results);
    
    %retrieve true and estimated spike times / window:
    [spiketimes_true_for_stats{si}, spiketimes_est_for_stats{si}, weights_true_for_stats{si}, weights_est_for_stats{si}, twin_for_stats(si, :)] = retrieve_spiketimes_for_stats(data, results);
    %note that only spiketimes within the relevant windows are returned
    
end

for si = 1:nseg
   
    s = segments_to_use(si);
    oedb.stats.bysegment{dataset}{neuron}(algind, s).totalspikes = sum(weights_true_for_stats{si});
    oedb.stats.bysegment{dataset}{neuron}(algind, s).totalT = diff(twin_for_stats(si, :));
    
end
oedb.stats.byneuron{dataset}(algind, neuron).totalspikes = sum(cellfun(@sum, weights_true_for_stats));
oedb.stats.byneuron{dataset}(algind, neuron).totalT = sum(twin_for_stats(:, 2) - twin_for_stats(:, 1));

%calculate correlation of spike counts
[~, ~, ~, csegment, cneuron] = oedb_masked_filter(spikecounts_true_for_stats, spikecounts_est_for_stats, oedb.opt.tol_sec, t_spikecounts_for_stats);
oedb.stats.byneuron{dataset}(algind, neuron).corr_spikecounts = cneuron;
for si = 1:nseg
    
    s = segments_to_use(si);
    oedb.stats.bysegment{dataset}{neuron}(algind, s).corr_spikecounts = csegment(si);
    
end

%calculate correlation of spike trains
[cneuron, ~, ~, ~, ~, ~, csegment] = point_process_corr_gfilt(spiketimes_true_for_stats, spiketimes_est_for_stats, oedb.opt.tol_sec, twin_for_stats, weights_true_for_stats, weights_est_for_stats);
oedb.stats.byneuron{dataset}(algind, neuron).corr_spiketimes = cneuron;
for si = 1:nseg
    
    s = segments_to_use(si);
    oedb.stats.bysegment{dataset}{neuron}(algind, s).corr_spiketimes = csegment(si);
    
end

%calculate detection and false positive rates from spike trains
[npairs, nspikes_true, nspikes_est, ~, ~, ntrueok_anymatch, nestok_anymatch] = test_spike_recon(spiketimes_true_for_stats, spiketimes_est_for_stats, oedb.opt.tol_sec, weights_true_for_stats, weights_est_for_stats);
oedb.stats.byneuron{dataset}(algind, neuron).det_spiketimes = sum(npairs) / sum(nspikes_true);
oedb.stats.byneuron{dataset}(algind, neuron).fp_spiketimes  = sum(nspikes_est - npairs) / sum(diff(twin_for_stats, 1, 2));
oedb.stats.byneuron{dataset}(algind, neuron).fp_spiketimes_anymatch  = sum(nspikes_est - nestok_anymatch) / sum(diff(twin_for_stats, 1, 2)); %FIXME do segments/dataset too
for si = 1:nseg
    
    s = segments_to_use(si);
    oedb.stats.bysegment{dataset}{neuron}(algind, s).det_spiketimes = npairs(si) / nspikes_true(si);
    oedb.stats.bysegment{dataset}{neuron}(algind, s).fp_spiketimes  = (nspikes_est(si) - npairs(si)) / diff(twin_for_stats(si, :));
    
end

%calculate detection and false positive rates from spike counts
[npairs, nspikes_true, nspikes_est, Tseg] = test_spike_recon_fromsc(spikecounts_est_for_stats, spikecounts_true_for_stats, t_spikecounts_for_stats, oedb.opt.tol_sec);
oedb.stats.byneuron{dataset}(algind, neuron).det_spikecounts = sum(npairs) / sum(nspikes_true);
oedb.stats.byneuron{dataset}(algind, neuron).fp_spikecounts  = sum(nspikes_est - npairs) / sum(Tseg);
for si = 1:nseg
    
    s = segments_to_use(si);
    oedb.stats.bysegment{dataset}{neuron}(algind, s).det_spikecounts = npairs(si) / nspikes_true(si);
    oedb.stats.bysegment{dataset}{neuron}(algind, s).fp_spikecounts  = (nspikes_est(si) - npairs(si)) / Tseg(si);
    
end

if all(cat(1, weights_true_for_stats{:}) == 1) && all(cat(1, weights_est_for_stats{:}) == 1)
    
    %analyze isolated single APs
    %FIXME do segments/dataset too
    
    [n_isolated_singles, ndetected] = test_singleAP_detection(spiketimes_true_for_stats, spiketimes_est_for_stats, oedb.opt.tol_sec);
    oedb.stats.byneuron{dataset}(algind, neuron).singleAP_rate_true = n_isolated_singles / sum(diff(twin_for_stats, 1, 2));
    oedb.stats.byneuron{dataset}(algind, neuron).singleAP_detrate = ndetected / n_isolated_singles;
    
    %analyze isolated bursts
    %FIXME do segments/dataset too
    
    [n_isolated_bursts, ndetected] = test_burst_detection(spiketimes_true_for_stats, spiketimes_est_for_stats, oedb.opt.tol_sec);
    oedb.stats.byneuron{dataset}(algind, neuron).burst_detrate = ndetected / n_isolated_bursts;

end


%calculate firing rates based on spike counts
oedb.stats.byneuron{dataset}(algind, neuron).fr_true_spikecounts = sum(cellfun(@sum, spikecounts_true_for_stats)) / sum(Tseg);
oedb.stats.byneuron{dataset}(algind, neuron).fr_est_spikecounts  = sum(cellfun(@sum, spikecounts_est_for_stats))  / sum(Tseg);
for si = 1:nseg
    
    s = segments_to_use(si);
    oedb.stats.bysegment{dataset}{neuron}(algind, s).fr_true_spikecounts = sum(spikecounts_true_for_stats{si}) / Tseg(si);
    oedb.stats.bysegment{dataset}{neuron}(algind, s).fr_est_spikecounts  = sum(spikecounts_est_for_stats{si})  / Tseg(si);
    
end

%calculate firing rates based on spike times
oedb.stats.byneuron{dataset}(algind, neuron).fr_true_spiketimes = sum(cellfun(@sum, weights_true_for_stats)) / sum(diff(twin_for_stats, 1, 2));
oedb.stats.byneuron{dataset}(algind, neuron).fr_est_spiketimes  = sum(cellfun(@sum, weights_est_for_stats))  / sum(diff(twin_for_stats, 1, 2));
for si = 1:nseg
    
    s = segments_to_use(si);
    oedb.stats.bysegment{dataset}{neuron}(algind, s).fr_true_spiketimes = sum(weights_true_for_stats{si}) / diff(twin_for_stats(si, :));
    oedb.stats.bysegment{dataset}{neuron}(algind, s).fr_est_spiketimes  = sum(weights_est_for_stats{si})  / diff(twin_for_stats(si, :));
    
end

%update dataset stats based on updated neuron stats
oedb = oedb_calc_dataset_stats(oedb, algind, dataset); %this might get called more than necessary but the computation is negligible so who cares?


function oedb = oedb_calc_dataset_stats(oedb, algind, datasetind)
fn = {'corr_spikecounts', 'corr_spiketimes', ...
    'det_spikecounts', 'det_spiketimes', ...
    'fp_spikecounts', 'fp_spiketimes', ...
    'fr_est_spikecounts', 'fr_est_spiketimes', 'fr_true_spikecounts', 'fr_true_spiketimes', ...
    'totalT', 'totalspikes' ...
    };

for k = 1:numel(fn)
    
    v = cat(2, oedb.stats.byneuron{datasetind}(algind, :).(fn{k}) );
    v(isnan(v)) = [];
    oedb.stats.bydataset(algind, datasetind).(fn{k})         = mean(v);
    oedb.stats.bydataset(algind, datasetind).([fn{k} '_sd']) =  std(v, 1); %don't do bias correction
    
end