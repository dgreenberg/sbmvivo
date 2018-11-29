function [r, alginfo] = apdet_MLspike(oerec, params, opts)
%[results, alginfo] = apdet_constrainedfoopsi(oerec, params, opts)
%
%https://github.com/MLspike/spikes/
assert(nargin == 3, '3 inputs required');

assert(exist('tps_mlspikes.m', 'file') == 2, 'tps_mlspikes.m was not found on the Matlab path. see https://github.com/MLspike/spikes/');
mlspikedir = fileparts(which('tps_mlspikes.m'));
brickdir = [mlspikedir filesep '..' filesep 'brick'];
assert(exist(brickdir, 'dir') == 7, 'the brick toolbox required for MLspike must be in the same folder as MLspike');

nrecs = numel(oerec);
ndatasegments = arrayfun(@(o) numel(o.data), oerec);
r = apdet_resultsstruct(ndatasegments);

for n = 1:nrecs    
    
    [f, it, dt,~,indicatorstring] = extract_data_fromoerec(oerec(n), 1:ndatasegments(n));
    
    assert(all(strcmpi(indicatorstring, indicatorstring{1})), 'not yet implemented for recordings of the same neuron with multiple fluorescence indicators');
    indicatorstring = indicatorstring{1};
    
    par = tps_mlspikes_wrapper('par');        
    par.drift.parameter = 0.01;
    par.display = 'none';  % don't pop up windows etc.
    par.dographsummary = false;
    
    %values from https://www.nature.com/articles/ncomms12190#supplementary-information
    switch lower(indicatorstring)
        case 'gcamp6s'
            
            par.a = 0.1130;
            par.tau = 1.87;
            par.pnonlin = [0.81 -0.056];
            
        case 'gcamp6f'
            
            par.a = 0.0341;
            par.tau = 0.76;
            par.pnonlin = [0.85 -0.006];
            
        case 'ogb'
            
            par.a = 0.049;
            par.tau = 0.78;
            
        otherwise
            error('unrecognized indicators: %s', indicatorstring);
    end
    
    par = repmat(par, [1 ndatasegments(n)]);
        
    f = cellfun(@(x) reshape(x, 1, []), f, 'uniformoutput', false);
    
    dt = nan(1, ndatasegments(n));
    for s = 1:ndatasegments(n)
        
        r{n}(s).spikecount_times = it{s} - dt(s) / 2;  % MLspike methods: "n_t (number of spikes between time t+1 and t)"
        par(s).dt = dt(s);
        r{n}(s).spiketimes_window = [it{s}(1) - dt(s) / 2, it{s}(end) + dt(s) / 2];  % it{1} since we called extract_data_fromoerec for this segment
        
    end
    
    %data for each dt must be processed separately
    alldt = unique(dt);
    spikecounts = cell(1, ndatasegments(n));
    for nextdt = alldt(:)'
        
        ii = dt == nextdt;
        spikecounts(ii) = tps_mlspikes_wrapper(f(ii), par(ii));
        
    end
    for s = 1:ndatasegments(n)
        
        spikecounts{s} = double(spikecounts{s});
        
    end

%     par = spk_autocalibration('par');
%     par.dt = dt;
%     [tauest aest sigmaest] = spk_autocalibration(f,par);
%     %ton .0702
%     %,'c0', 0.36 ...
%     %    saturation 7.9e-3
%     %hill exponent 1.84

    for s = 1:ndatasegments(n)
        
        r{n}(s).opts = opts;
        r{n}(s).params = par(s);
        
        r{n}(s).spikecounts = spikecounts{s};
        r{n}(s).outputvars = struct();
        
        spiketimes = [];
        for k = find(spikecounts{s}(:)' ~= 0)
            
            tmidpoint = r{n}(s).spikecount_times(k); %midpoint between two scan times
            if spikecounts{s}(k) == 1
                
                spiketimes(end + 1, 1) = tmidpoint; %#ok<AGROW>
                
            else
                
                spiketimes = [spiketimes; tmidpoint - dt(s) / 2 + dt(s) * (0.5 + (0:spikecounts{s}(k) - 1)') / spikecounts{s}(k)]; %#ok<AGROW>
                
            end
            
        end
        r{n}(s).spiketimes = spiketimes;
        
    end
    
end

alginfo = empty_apdetalginfo;
alginfo.name = 'MLspike';
alginfo.shortname = 'MLspike';
alginfo.ismaximumlikelihood = true;
alginfo.isdeterministic = true;
alginfo.algorithm_version = 0.01;