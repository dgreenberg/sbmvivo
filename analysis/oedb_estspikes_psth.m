function [psth, bc, nspikes, meantdiff, tdiff_w, bias] = oedb_estspikes_psth(oedb, fileinfo, dataset, algind, hist_binsize, hist_win, mindt, maxdt)
if ~exist('hist_win', 'var') || isempty(hist_win)

    hist_win = [-0.5 0.5];
    
end
if ~exist('hist_binsize', 'var') || isempty(hist_binsize)

    hist_binsize = 0.05;
    
end
if ~exist('mindt','var') || isempty(mindt), mindt = -inf; end
if ~exist('maxdt','var') || isempty(maxdt), maxdt = inf; end
isolation_win = [min(-1, hist_win(1)) max(0.5, hist_win(2))];

oerec = fetch_dataset_oerecarray(oedb, fileinfo, dataset);

nneurons = numel(oerec);
nspikes = zeros(1, nneurons);
tdiff_w = cell(1, nneurons);
for n = 1:nneurons
    
    dt_eachseg = arrayfun(@(v) diff(v.t(1:2)) * v.t_scale_factor, oerec(n).data);
    
    for d = 1:numel(oerec(n).data)
        
        if dt_eachseg(d) < mindt || dt_eachseg(d) > maxdt
            
            continue;
            
        end
        
        st_true = oerec(n).data(d).spiketimes;
        results = fetch_results(oedb, fileinfo, dataset, n, d, algind);
        if ~any(isnan(results.spiketimes_window))
            
            st_est = results.spiketimes(:);
            w_est = ones(size(st_est));
            twin = results.spiketimes_window;
            
        elseif ~any(isnan(results.spikecounts)) && ~isempty(results.spikecounts)
            
            [st_est, w_est, twin] = convert_spikecounts_to_weighted_spiketrain(results.spikecounts, results.spikecount_times);
            
        else
            
            continue;
            
        end
        st_true = st_true(st_true > twin(1) & st_true < twin(2));
        
        for j = 1:numel(st_true)
            
            nextst = st_true(j);
            if any(st_true < nextst & st_true > nextst + isolation_win(1)) || any(st_true > nextst & st_true < nextst + isolation_win(2))
                
                continue;
                
            end
            
            ii = st_est > nextst + hist_win(1) & st_est < nextst + hist_win(2);
            tdiff_w{n} = [tdiff_w{n}; st_est(ii) - nextst, w_est(ii)];
            nspikes(n) = nspikes(n) + 1;  % counts true spikes
            
        end
        
    end
    
end

ed = (floor(hist_win(1) / hist_binsize):ceil(hist_win(2) / hist_binsize)) * hist_binsize;
bc = ed(1:end - 1) + diff(ed(1:2)) / 2;
psth = zeros(numel(bc), nneurons);
[bias, meantdiff] = deal(nan(1, nneurons));
for n = 1:nneurons
    
    if ~nspikes(n), continue; end
    [~, jj] = histc(tdiff_w{n}(:, 1), ed);
    for u = 1:numel(bc)
        
        psth(u, n) = sum(tdiff_w{n}(jj == u, 2));
        
    end
    
    meantdiff(n) = sum(abs(tdiff_w{n}(:, 1)) .* tdiff_w{n}(:, 2)) / sum(tdiff_w{n}(:, 2));
    bias(n) = sum(tdiff_w{n}(:, 1) .* tdiff_w{n}(:, 2)) / sum(tdiff_w{n}(:, 2));
    
end
psth = bsxfun(@rdivide, psth, nspikes);