function [Pout, opts, firingrates] = globalfit(varargin)
% [Pout, opts, firingrates] = ...
%     globalfit(f, it, st, neuronindex, neuronnames, basefigsavedir, P0, ...
%     sessionindex, fixed_params_fromuser, focalplaneindex, modelname, target_substepsize, maxtotalt_hours, ...
%     savefigs, varsfile, opts)

%fixme convert some inputs to opts fields
%FIXME should estimate fluorescence noise, will require nA2D

%parse inputs
[f, it, st, neuronindex, neuronnames, basefigsavedir, P0, sessionindex, fixed_params_fromuser, focalplaneindex, ...
    modelname, target_substepsize, maxtotalt_hours, savefigs, varsfile, opts, savevars, neuronnames_eachsegment] = parseinputs_fitautoR(varargin{:});

forig = f;
opts = sbm.fit.assignopts(opts);
nseg = numel(f);

[simulation_function, default_parameters, parameter_names, statenames, always_fixed_params, params_to_print] = initialize_fit(opts);
fixed_params = union(always_fixed_params, fixed_params_fromuser);
parameter_names = setdiff(parameter_names, fixed_params);
params_to_print = setdiff(params_to_print, fixed_params);

P0 = apply_default_params_autoR(P0, default_parameters, nseg, opts);

%calculate baseline fluorescence using AP times
[fBL, blmask] = get_baseline_forfit(f, it, st, opts.tau_s0, opts.verbose);
%normalize f and fBL by the median fBL value for each session
medf_orig = nan(1, numel(sessionindex));
for si = reshape(unique(sessionindex), 1, [])
    
    m = abs(median(cat(1, fBL{sessionindex == si}))); %if fBL has mostly negative values, don't flip its sign!!
    medf_orig(sessionindex == si) = m;
    for ii = reshape(find(sessionindex == si), 1, [])
        
        f{ii} = f{ii} / m;
        fBL{ii} = fBL{ii} / m;
        
    end
    
end

%split data segments into subsegments, throwing out data far from all spikes to reduce the amount of ODE solving we have to do
[fsplit, itsplit, stsplit, fBLsplit, timeind, segmentind_split, blmasksplit] = split_segs( ...
    f, it, st, fBL, blmask, opts.perispikewin);
nsubseg = numel(fsplit);
segment_used = ismember(1:nseg, segmentind_split);
nneurons_used = numel(unique(neuronindex(segmentind_split)));
assert(any(segment_used), 'no usable data');

%renumber indices so they cover consecutive integers starting with 1
neuronindex_orig = neuronindex; focalplaneindex_orig = focalplaneindex; %save these so we can create parameter outputs of the right sizes while including unused data
[neuronindex, sessionindex, focalplaneindex, nneurons, nsessions, nfocalplanes] = renumber_indices(neuronindex, sessionindex, focalplaneindex, segment_used);

%determine step size, round each spike to the nearest modeling timestep, and get the time course for states/spikes
[dt_eachsubseg, n_extrabins_eachsubseg, bin_edges, t_states, stepsize_eachsubseg, nsteps_eachsubseg] = timinginfo( ...
    itsplit, target_substepsize, nsubseg, stsplit);

[states, states_eq] = initialize_states(P0, t_states, opts.modelname);

if opts.verbose > 5
    
    fprintf('Using model %s.\n', opts.modelname);
    fprintf('Maximum simulation step size is %g sec.\n', target_substepsize);
    
end

%build objective function and get info on param vector elements, and initial    paramter vector
[p0, pi, p_seglist, pmin, pmax, isobs] = parameterization(opts, P0, parameter_names, nseg, neuronindex, sessionindex, focalplaneindex);
parameter_names(~ismember(parameter_names, fieldnames(pi))) = [];

[nspikes_eachneuron, T_eachneuron] = deal(zeros(1, max(neuronindex_orig)));
for si = 1:nseg %only keep data with imaging, as we dont' want to include cell dying etc. and spikes might not have been carefully checked
    
    n = neuronindex_orig(si);
    dt = diff(it{si}(1:2));
    nspikes_eachneuron(n) =  nspikes_eachneuron(n) + sum(st{si} >= it{si}(1) - dt & st{si} <= it{si}(end));
    T_eachneuron(n) = T_eachneuron(n) + it{si}(end) - it{si}(1) + dt;
    
end
firingrates = nspikes_eachneuron ./ T_eachneuron;

spikecounts = getspikecounts(stsplit, P0(1).spiketshift, bin_edges);

%the Pfunc gives us a paramter structure array from a parameter vector
Pfunc = @(p) pvec2Pstruct(P0, parameter_names, pi, p, pmin, opts);

ofunc_split = @(p) error_function( ...
    Pfunc, simulation_function, ...
    fsplit, fBLsplit, spikecounts, stepsize_eachsubseg, nsteps_eachsubseg, ...
    n_extrabins_eachsubseg, p, pmin, isobs, opts.symgrad, ...
    segmentind_split, p_seglist, opts.gradreldiff, opts.gradmindiff, neuronindex(segmentind_split), focalplaneindex(segmentind_split), opts, ...
    states, states_eq, true, 1:nsubseg);

commandwindow; %so we can see the program output while running the fitting procedure
pause(0.2);

%test function at initial parameters, and time it
tic;
[e0, g0, fhatsplit0, statessplit0, esegsplit0, R0, states_eq_eachsubseg0, eneuron0] = ofunc_split(p0);  % this currently causes a gradient evaluation. imperfect, but the line search seems to as well
tfg = toc;

if opts.verbose > 5
    
    fprintf('Function + gradient evaluation time at starting parameters: %g sec.\n', tfg);
    
end

maxiter = min(10000, ceil(maxtotalt_hours * 3600 / tfg));
maxfunevals = maxiter * 20;  %allow up to 20 function evaluations per iteration

if opts.verbose > 5
    
    fprintf('Maximum number of iterations: %d\n', maxiter);
    
end

%optimize
displaymode = 'iter';
if opts.verbose < 7, displaymode = 'off'; end
figh_opf = figure('visible','off');  %store data in invisible figure, surely there's a better way to do this
opf = @(x, optimValues, state) canonlinfit_opt_outputfcn(x, optimValues, state, figh_opf);
n1 = now;
if opts.uselog
    
    if any(pmax ~= inf)
        
        warning('in log mode, maximum values for all parameters are ignored');
        
    end
    
    oopts = optimset('maxfunevals',maxfunevals,'maxiter',maxiter,'tolfun',opts.tolfun,'tolx',opts.tolx,'largescale','off','display',displaymode,'gradobj','on','OutputFcn',opf);
    
    [p,e,ef,op,grad] = fminunc(ofunc_split, p0, oopts);
    
else
    
    oopts = optimset('maxfunevals',maxfunevals,'maxiter',maxiter,'tolfun',opts.tolfun,'tolx',opts.tolx,'display',displaymode,'gradobj','on','OutputFcn',opf,...
        'algorithm','sqp','tolcon',1e-10, 'MaxSQPIter', 300);
    
    [p,e,ef,op,~,grad] = fmincon(ofunc_split, p0, [], [], [], [], pmin, pmax, [], oopts);
    
end
n2 = datevec(now - n1);
if opts.verbose > 5
    fprintf('total time : %f sec.\n',n2(end-2:end) * [3600 60 1]')
end
allp = getappdata(figh_opf, 'allp');
alle = getappdata(figh_opf, 'alle');
close(figh_opf);

%get state dynamics etc. after convergence
[e, grad, fhatsplit, statessplit, esegsplit, R, states_eq_eachsubseg, eneuron] = ofunc_split(p);
%at this point statessplit does yet match statenames

%augment states and states_eq with calcium-free states of indicator and of each buffer
P = Pfunc(p);
P_eachsubseg = P(segmentind_split);
[statessplit, states_eq_eachsubseg] = augment_states(statessplit, states_eq_eachsubseg, P_eachsubseg, opts.modelname);
%statessplit now matches statenames

%print parameter changes
if opts.verbose > 5
    canonlinfit_print_changes(e0 / nneurons_used, e / nneurons_used, Pfunc, p, p0, params_to_print, pi, R0, R);
end

%put results back together. FIXME: combine multiple segs in a session too (would have to take another input to do that)
[fBL_combo, blmask_combo, fhat_combo, states_combo, t_states_combo, eseg_combo] = combine_segs( ...
    f, fBLsplit, blmasksplit, fhatsplit, statessplit, segmentind_split, timeind, t_states, esegsplit);
[Pout, L] = Pfunc(p);

% compute RMS errors
ntimepoints_eachseg_combo = cellfun(@numel, f);
rms_combo = sqrt(eseg_combo(:) ./ ntimepoints_eachseg_combo(:));

[Pout, Rind] = canonlinfit_update_Pout(Pout, focalplaneindex, neuronindex, focalplaneindex_orig, neuronindex_orig, opts, R);

if savefigs
    
    canonlinfit_save_figures(Pout, segmentind_split, neuronnames_eachsegment, states_combo, states_eq_eachsubseg, t_states_combo, ...
        it, fBL, f, fhat_combo, blmask_combo, rms_combo, R, Rind, L, pi, statenames, opts.perispikewin, basefigsavedir);
    
end
if savevars
    
    clear esegfig;  % don't save a function handle
    
    clear states states_combo states statessplit statessplit0 bin_edges t_states spikecounts t_states_combo varargin;  % don't save huge vars, we can regenerate them if needed from the pfunc/pval
    save([basefigsavedir varsfile]);
    
end


function [n, s, f, nn, ns, nf] = renumber_indices(n0, s0, f0, segment_used)
[~, ~, n_used] = unique(n0(segment_used));
[~, ~, s_used] = unique([n0(segment_used) s0(segment_used)], 'rows');
[~, ~, f_used] = unique([n0(segment_used) f0(segment_used)], 'rows');

[n, s, f] = deal(zeros(size(n0)));
n(segment_used) = n_used;
s(segment_used) = s_used;
f(segment_used) = f_used;

nn = max(n);
ns = max(s);
nf = max(f);


function [f, it, st, neuronindex, neuronnames, basefigsavedir, P0, sessionindex, fixed_params_fromuser, focalplaneindex, ...
    modelname, target_substepsize, maxtotalt_hours, savefigs, varsfile, opts, savevars, neuronnames_eachsegment] = parseinputs_fitautoR( ...
    f, it, st, neuronindex, neuronnames, basefigsavedir, P0, sessionindex, fixed_params_fromuser, focalplaneindex, ...
    modelname, target_substepsize, maxtotalt_hours, savefigs, varsfile, opts)
if nargin < 3
    
    error('f, it and st must be provided');
    
end
if isnumeric(f)
    f = {f}; st = {st}; it = {it};
end
nseg = numel(f);
for ii = 1:nseg, f{ii} = reshape(f{ii}, [], 1); end

if ~exist('maxtotalt_hours', 'var') || isempty(maxtotalt_hours)
    
    maxtotalt_hours = 20;
    
end
if ~exist('target_substepsize', 'var') || isempty(target_substepsize)
    
    target_substepsize = 1e-3;
    
end
if exist('basefigsavedir','var') && ~isempty(basefigsavedir) && exist(basefigsavedir,'dir') == 7
    
    if ~exist('savefigs', 'var') || isempty(savefigs)
        
        savefigs = false;
        
    end
    if ~exist('varsfile', 'var') || isempty(varsfile)
        
        varsfile = '';
        
    end
    savevars = exist('varsfile', 'var') && ~isempty(varsfile);
    
else
    
    [savefigs, savevars] = deal(false);
    
end

if ~exist('sessionindex','var') || isempty(sessionindex)
    
    sessionindex = 1:numel(f);
    
else
    
    assert(numel(sessionindex) == numel(f));
    
end
if ~exist('neuronindex','var')
    
    neuronindex = ones(nseg, 1);
    
else
    
    assert(numel(neuronindex) == numel(f));
    
end
if ~exist('focalplaneindex','var') || isempty(focalplaneindex)
    
    focalplaneindex = 1:numel(f);
    
else
    
    assert(numel(focalplaneindex) == numel(f));
    
end
if ~exist('neuronnames','var') || isempty(neuronnames), neuronnames = repmat({''}, 1, max(neuronindex)); end
if ~exist('P0','var'), P0 = []; end
if ~exist('opts', 'var') || isempty(opts), opts = struct(); end
if ~exist('fixed_params_fromuser', 'var') || isempty(fixed_params_fromuser)
    
    fixed_params_fromuser = {};
    
else
    
    fixed_params_fromuser = fixed_params_fromuser(:)';
    
    if ismember('k50', fixed_params_fromuser)
        
        fixed_params_fromuser = [fixed_params_fromuser {'k50_first' 'dk50'}];
        if ismember('tau_mm', fixed_params_fromuser)
            
            fixed_params_fromuser = [fixed_params_fromuser {'koff' 'koff'}];
            
        end
        
    end
    
end
neuronnames_eachsegment = neuronnames(neuronindex); %so that we save the figure files with the correct names


function P0 = apply_default_params_autoR(P0, default_parameters, nseg, opts)
if isempty(P0)
    
    P0 = repmat(default_parameters, 1, nseg);
    
else
    
    if numel(P0) == 1
        
        P0 = repmat(P0, 1, nseg);
        
    else
        
        assert(numel(P0) == nseg, 'P must have at the same number of elements as there are data segments');
        
    end
    if isfield(P0,'k50')
        for u = 1:nseg
            
            P0(u).k50_first = P0(u).k50(1);
            if numel(P0(u).k50) > 1
                
                P0(u).dk50 = diff(P0(u).k50);
                
            end
            
        end
    end
    
    fn = fieldnames(default_parameters);
    for u = 1:nseg
        
        for k = 1:numel(fn)
            
            if ~isfield(P0, fn{k}) || isempty(P0(u).(fn{k}))
                
                P0(u).(fn{k}) = default_parameters.(fn{k});
                
            end
            
        end
        
    end
    
end
P0 = update_dependent_params_autoR(P0, opts, 1:nseg);