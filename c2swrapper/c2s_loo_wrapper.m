function [ns, r] = c2s_loo_wrapper(oerec)
inputfilename = [fileparts(mfilename('fullpath')) filesep 'c2s_data_file_tmp.mat'];
outputfilename = [fileparts(mfilename('fullpath')) filesep 'c2s_loo_results_tmp.mat'];
oerec2c2sfile(inputfilename, oerec);

commandstr = sprintf('c2s leave-one-out %s %s', inputfilename, outputfilename);

assert(system(commandstr) == 0, 'command failed');
r = load(outputfilename, '-mat');

nneurons = numel(oerec);
nseg = arrayfun(@(v) numel(v.data), oerec);
ns = cell(1, nneurons);
ii = 0;
for n = 1:nneurons    
    for s = 1:nseg(n)
        ii = ii + 1;
        t = double(oerec(n).data(s).t) * oerec(n).data(s).t_scale_factor + oerec(n).data(s).t_offset;
        dt = median(diff(t));
        
        %c2s uses scipy.signal.resample to resample
        %https://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.signal.resample.html        
        first_prediction_time = t(1);
        nf = numel(t);
        npred = numel(r.data{ii}.predictions);        
        dx_fac = nf / npred;
        last_prediction_time_in_origsamples = (npred - 1) * dx_fac;
        prev_sample = find((1:nf) <= last_prediction_time_in_origsamples, 1, 'last');
        assert(numel(prev_sample) == 1, 'failed to align timing');        
        extra_time_in_samples = last_prediction_time_in_origsamples - prev_sample;
        last_prediction_time = t(prev_sample) + extra_time_in_samples * dt;
        t_predictions = linspace(first_prediction_time, last_prediction_time, npred);        
        
        histogram_bin_edges = [t(1) - dt; reshape(t, [], 1)];
        [~, bin_index] = histc(t_predictions, histogram_bin_edges);
        
        ns{ii} = nan(numel(t), 1);
        for k = 1:numel(t)
            ns{ii}(k) = sum(r.data{ii}.predictions(bin_index == k));
        end
        ns{ii}(end) = ns{ii}(end) + sum(r.data{ii}.predictions(bin_index == numel(t) + 1)); %last histogram bin contains samples that precisely match the last histogram bin edge        
    end
end

delete(inputfilename);
delete(outputfilename);