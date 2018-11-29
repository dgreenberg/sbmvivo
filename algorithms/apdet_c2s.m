function [r, alginfo] = apdet_c2s(oerec, params, opts)
%[results, alginfo] = apdet_c2s(oerec, params, opts)
assert(nargin == 3, '3 inputs required');

if isempty(oerec)  % the case when querying algorithm info
    
    [commandresult, ~] = system('c2s -h');
    assert(commandresult == 0, 'c2s is not installed on this system');
    
end

nrecs = numel(oerec);
ndatasegments = arrayfun(@(o) numel(o.data), oerec);
r = apdet_resultsstruct(ndatasegments);

modelsdir = [fileparts(which('c2s_predict_wrapper')) filesep 'models' filesep];

for n = 1:nrecs
    
    assert(all(arrayfun(@(d) all(d.f_mask), oerec(n).data)), 'not implemented for data with gaps');
    nextparams = params{n}{1};
    
    [~,it,dt,~,indicatorstring,~] = extract_data_fromoerec(oerec(n), 1:ndatasegments(n));
    
    if ~isfield(nextparams, 'modelfile') || isempty(nextparams.modelfile)
        
        %fixme let default model depend on indicatorstring
        nextparams.modelfile = ''; %default model
        modelfile_fullpath = '';
        %nextparams.modelfile = [modelsdir 'gcamp6sv1mouse.xpck']; %fixme deal with multiple indicators etc.
        
    else
        
        modelfile_fullpath = [modelsdir filesep nextparams.modelfile];
        assert(exist(modelfile_fullpath, 'file') == 2, 'model file not found in %s', modelsdir);
        
        for s = 2:ndatasegments(n)
            
            assert(strcmpi(nextparams.modelfile, params{n}{s}.modelfile), 'same model must be used for all a neuron''s segments');
            
        end
        
    end
    
    %run c2s on all the neuron's data at once:
    rneuron = c2s_predict_wrapper(oerec(n), modelfile_fullpath);
    %note that because we have called c2s predict on one neuron at a time, we don't have to worry about cell_num etc.
    
    %format outputs
    for s = 1:ndatasegments(n)
        
        r{n}(s).params = nextparams;
        r{n}(s).opts = opts;
        
        r{n}(s).spikecounts = rneuron{s}.predictions(:);
        
        %c2s uses scipy.signal.resample to resample
        %https://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.signal.resample.html        
        first_prediction_time = it{s}(1);
        nf = numel(it{s});
        npred = numel(rneuron{s}.predictions);        
        dx_fac = nf / npred;
        last_prediction_time_in_origsamples = (npred - 1) * dx_fac;
        prev_sample = find((1:nf) <= last_prediction_time_in_origsamples, 1, 'last');
        assert(numel(prev_sample) == 1, 'failed to align timing');        
        extra_time_in_samples = last_prediction_time_in_origsamples - prev_sample;
        last_prediction_time = it{s}(prev_sample) + extra_time_in_samples * dt(s);
        r{n}(s).spikecount_times = linspace(first_prediction_time, last_prediction_time, npred)'; 
        
    end
    
end

alginfo = empty_apdetalginfo;
alginfo.name = 'c2s';
alginfo.shortname = 'c2s';
alginfo.ismaximumlikelihood = false;
alginfo.isdeterministic = true;
%alginfo.algorithm_version = ; fixme