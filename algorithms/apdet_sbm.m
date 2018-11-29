function [r, alginfo] = apdet_sbm(oerec, params, opts)
%[results, alginfo] = apdet_sbm(oerec, params, opts)
assert(nargin == 3, '3 inputs required');

nrecs = numel(oerec);
ndatasegments = arrayfun(@(o) numel(o.data), oerec);
r = apdet_resultsstruct(ndatasegments);

for n = 1:nrecs
    
    [F,it,dt,~,indicatorstring,nA2D] = extract_data_fromoerec(oerec(n), 1:ndatasegments(n));
    
    [opts.settingsindex, opts.sessionindex, opts.darkoffsetcorrected, opts.fdcindex, opts.focalplaneindex] = oerec_datagroupindices_singleneuron(oerec(n)); %FIXME remove these from opts
    
    assert(all(strcmpi(indicatorstring, indicatorstring{1})), 'sbm is not yet implemented for recordings of the same neuron with multiple fluorescence indicators');
    indicatorstring = indicatorstring{1};
    
    assert(all(arrayfun(@(d) all(d.f_mask), oerec(n).data)), 'sbm is not yet implemented for data with gaps');
    
    P = cat(2, params{n}{:});
    
    opts.darkoffsetcorrected = true; %hack, fixme
    
    %run sbm-based inference on all the neuron's data at once:
    [M, P, opts, lik, prevP, stimx, moments, prevlik] = sbm.run_alg(F, dt, indicatorstring, opts, P, [], [], nA2D, []);
    
    %format outputs
    for s = 1:ndatasegments(n)
                
        r{n}(s).params = P(s);
        r{n}(s).opts = opts;
        
        r{n}(s).spikecounts = moments(s).n_mean;
        %times for spikecounts are shifted back one substep compared to the times for binding states etc.
        r{n}(s).spikecount_times = reshape(bsxfun(@plus, dt(s) * (0:moments(s).nsteps - 1)' / moments(s).nsteps, reshape(it{s},1,[]) - dt(s)),[],1);
        
        if isfield(M, 'gst') && ~any(isnan(M(s).gst.gfit))
            
            nframes = numel(it{s});
            nspikes = numel(M(s).gst.mu);
            st = nan(nspikes, 1);
            
            %the states and observations have to be shifted back by stepsize to get spike times:
            stepsize = dt(s) / moments(s).nsteps;
            st_steps = M(s).gst.mu / stepsize - 1;  % get the times of non-discretized detected APs in simulation steps. zero is exactly one frame (dt) before the first fluorescence measurement
            st_frames = st_steps / moments(s).nsteps;  % get the times of non-discretized detected APs in frames. zero is exactly one frame (dt) before the first fluorescence measurement
            st_frames_int = floor(st_frames);
            st_frames_frac = st_frames - st_frames_int;
            
            before_first_measurement = st_frames < 1;
            st(before_first_measurement) = it{s}(1) + (st_frames(before_first_measurement) - 1) * dt(s);
            during_measurement = st_frames >= 1 & st_frames <= nframes;
            st(during_measurement) = it{s}(st_frames_int(during_measurement)) + st_frames_frac(during_measurement) * dt(s);
            
            after_last_measurement = st_frames > nframes;
            st(after_last_measurement) = it{s}(end) + (st_frames(after_last_measurement) - nframes) * dt(s);
            
            r{n}(s).spiketimes = st;
            %extend window half a step before/after first/last possible time for a simulated spike:
            r{n}(s).spiketimes_window = [it{s}(1) - dt(s), it{s}(end) - stepsize] + 0.5 * stepsize * [-1 1];
            
        end
        
        r{n}(s).outputvars = orderfields(struct( ...
            'M', M(s), ...
            'prevP', prevP(:, s), ...
            'lik', lik(s), ...
            'prevlik', prevlik(:, s), ...
            'moments', moments(s) ...
            ));
        
    end
    
end

alginfo = empty_apdetalginfo;
alginfo.name = 'Sequential Binding Model';
alginfo.shortname = 'sbm';
alginfo.ismaximumlikelihood = false;
alginfo.isdeterministic = false;
alginfo.algorithm_version = sbm.version;