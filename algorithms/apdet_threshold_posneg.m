function [r, alginfo] = apdet_threshold_posneg(oerec, params, opts)
%[results, alginfo] = apdet_threshold_posneg(oerec, params, opts)
assert(nargin == 3, '3 inputs required');

nrecs = numel(oerec);
ndatasegments = arrayfun(@(o) numel(o.data), oerec);
r = apdet_resultsstruct(ndatasegments);

default_params = orderfields(struct( ...
    'target_false_positive_rate', 0.05, ...
    'baseline_percentile', 0.08, ...
    'baseline_window', 15, ... seconds. for calculating baseline
    'sigma_est_window', 1.0, ... seconds. we calculate sigma once for each window, then take the median
    'sigma_thresh_start', 4, ...
    'sigma_thresh_end', 2 ...
    ));
default_fn = fieldnames(default_params);

for n = 1:nrecs
    
    [dff, local_sigma_vals, fBL] = deal(cell(1, ndatasegments(n)));
    
    [fraw, it] = extract_data_fromoerec(oerec);
    medfBL = nan(1, ndatasegments(n));
    
    for s = 1:ndatasegments(n)
        
        r{n}(s).opts = opts;
        r{n}(s).params = params{n}{s};
        
        for k = 1:numel(default_fn)
            if ~isfield(r{n}(s).params, default_fn{k})
                
                r{n}(s).params.(default_fn{k}) = default_params.(default_fn{k});
                
            end
        end        
        
        dt = diff(it{s}(1:2));
        fraw = reshape(fraw, [], 1);
        
        fBL{s} = percentile_filter_windowed(fraw{s}, ceil(r{n}(s).params.baseline_window / dt), r{n}(s).params.baseline_percentile);
        medfBL(s) = median(fBL{s});
        
        dff{s} = (fraw{s} - fBL{s}) ./ medfBL(s);  % in addition to subtracting baseline, we normalize so that transient height will be the same across segments
        
        sigma_winsize_bins = round(r{n}(s).params.sigma_est_window / dt);
        assert(sigma_winsize_bins >= 5, 'temporal resolution too low to for s.d. window size of %f', r{n}(s).params.sigma_est_window);
        nwins = floor(numel(dff{s}) / sigma_winsize_bins);
        local_sigma_vals{s} = std(reshape(dff{s}(1:nwins * sigma_winsize_bins), nwins, []));
        
    end
    
    all_sigma_vals = cat(2, local_sigma_vals{:});
    assert(numel(all_sigma_vals) >= 10, 'not enough data to estimate sigma');
    sigma = median(all_sigma_vals);  % in dff
    
    thresh_start = r{n}(s).params.sigma_thresh_start * sigma; % in dff
    thresh_end = r{n}(s).params.sigma_thresh_end * sigma; % in dff
    
    for s = 1:ndatasegments(n)
        
        r{n}(s).spikecount_times = it{s};
        r{n}(s).spikecounts = zeros(numel(it{s}), 1);
        r{n}(s).outputvars = orderfields(struct('dff', dff{s}, 'sigma', sigma, 'fBL', fBL{s}));
        
        offset = 0;
        while true
            
            a = find(dff{s}(offset + 1:end) >= thresh_start, 1) + offset;
            if isempty(a), break; end
            offset = a;
            b = find(dff{s}(offset + 1:end) < thresh_end, 1) + offset;
            if isempty(b), b = numel(dff{s}); end
            offset = b;
            [~, j] = max(dff{s}(a:b));
            r{n}(s).spikecounts(a + j - 1) = 1;                       
            
        end
        
        r{n}(s).spiketimes = it{s}(r{n}(s).spikecounts > 0);
        dt = diff(it{s}(1:2));
        r{n}(s).spiketimes_window = [it{s}(1) - dt / 2, it{s}(end) + dt / 2];
        
    end
    
end

alginfo = empty_apdetalginfo;
alginfo.name = 'Threshold for positive-going events';
alginfo.shortname = 'pospeakthresh';
alginfo.ismaximumlikelihood = true;
alginfo.isdeterministic = true;
alginfo.algorithm_version = 0.01;
