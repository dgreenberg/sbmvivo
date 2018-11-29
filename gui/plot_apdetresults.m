function plot_apdetresults(results, dt, f_axes, spike_axes, algname, state_axes, normalize_states)
if ~any(isnan(results.spiketimes_window)) %plot spike times
    
    plotspiketrain(results.spiketimes, spike_axes, 0, 0.9);

end
if ~isempty(results.spikecount_times)
    
    line(results.spikecount_times, results.spikecounts, 'color', 'b', 'parent', spike_axes);
    
end
if ~exist('normalize_states', 'var') || isempty(normalize_states)
    
    normalize_states = false;
    
end

switch lower(algname)
    case 'sbm'
        
        if ~isempty(fieldnames(results.outputvars)) && ~isempty(state_axes)
            
            show_sbm_graphics(f_axes, spike_axes, state_axes, results.outputvars, dt, results.spikecount_times, results.params, results.opts, normalize_states, results.spiketimes);
            
        end
        
end

if ~isempty(results.spikecounts) || ~isempty(results.spiketimes)
    
    set(spike_axes,'ylimmode','auto');
    YL = get(spike_axes,'ylim');
    set(spike_axes,'ylim',[-0.1 max(YL(2), 3)]);    
    
else
    
    YL = get(spike_axes,'ylim');
    set(spike_axes,'ylim',[-0.1 max(YL(2), 3)]);
    
end


function show_sbm_graphics(f_axes, spike_axes, state_axes, outputvars, dt, estns_t, params, opts, normalize_states, spiketimes)
if ~isfield(outputvars, 'M') || ~isfield(outputvars, 'moments'), return; end
if ~isfield(params, 'S'), return; end
M = outputvars.M; moments = outputvars.moments;

%we shift back it_sub for spike moments, but not for states
t_states = estns_t + dt / moments.nsteps;
t_observations_used = t_states(moments.nsteps:moments.nsteps:end); %only observations used if we used the option "only process onscreen data"

if isfield(params, 'fdc')
    fdc = params.fdc;
else
    fdc = 0;
end

Fmean = moments.gp_mean' + fdc;
Fsd = moments.Fpred_sd'; %total posterior variance of F given hidden states, due to both observation noise and uncertainty in state values

if strcmpi(opts.model, '3s')
    
    [s0, x0, ~] = sbm.model.equilibriumstates(params, opts);%calculate equilibrium states
    brightness_eq = s0 + params.alpha_i * x0 + params.R;
    
elseif strcmpi(opts.model, '5s2b')
    
    [seq, ~] = sbm.model.equilibriumstates(params, opts);
    brightness_eq = (params.dbrightness(:)' * seq / params.S + 1) * (1 + params.R);
    
else
    
    error('unrecognized model: %s', opts.model);
    
end

FBLmean = moments.expr_mean' * brightness_eq + fdc;
FBLlow = FBLmean - moments.expr_sd' * brightness_eq;
FBLhigh = FBLmean + moments.expr_sd' * brightness_eq;

nsmean = moments.n_mean';
nslow = max(0, nsmean - moments.n_sd');
nshigh = nsmean + moments.n_sd';

Fintpatch = patch([t_observations_used; flipud(t_observations_used); t_observations_used(1)], [Fmean - Fsd; flipud(Fmean + Fsd); Fmean(1) - Fsd(1)],'c','parent',f_axes,'edgecolor','none');
Fmeanline = line(t_observations_used, Fmean, 'parent', f_axes,'color','b','linewidth',2);

patchx_aps = reshape([estns_t; flipud(estns_t); estns_t(1)],[],1); %x vals for states etc.
nsintpatch = patch(patchx_aps, [nslow; flipud(nshigh); nslow(1)],'c','parent',spike_axes,'edgecolor','none');

ch = get(spike_axes,'children');
set(spike_axes,'children',[reshape(setdiff(ch, reshape(nsintpatch,[],1)),[],1); nsintpatch]);

if isfield(opts,'gfitst') && opts.gfitst && isa(M, 'struct') && isfield(M,'gst')
    
    gwin = 0.4;
    
    stepsize = diff(M.gst.tt(1:2));
    gwin_bins = ceil(gwin / (stepsize/10));
    
    dx = (-gwin_bins:gwin_bins)' * (stepsize/10);
    x = bsxfun(@plus, dx, spiketimes(:)');
    dxnorm = bsxfun(@rdivide, dx, reshape(M.gst.sd, 1, []));
    
    y = bsxfun(@rdivide, exp(-0.5 * dxnorm .^ 2), sqrt(2 * pi) * reshape(M.gst.sd, 1, [])) * (sqrt(2 * pi) * stepsize);
    y(abs(dxnorm) > 3) = NaN; %only show 3 s.d.
    
    line(estns_t, M.gst.gfit, 'parent',spike_axes,'color', [0 .6 0],'linewidth',2);
    line(reshape([x; nan(1, size(x, 2))], [], 1), reshape([y; nan(1, size(y, 2))], [], 1),'color','m','parent',spike_axes);
            
end

fblpatch = patch([t_observations_used; flipud(t_observations_used); t_observations_used(1)], [FBLlow; flipud(FBLhigh); FBLlow(1)],[1 .6 .6],'parent',f_axes,'edgecolor','none');
fblline = line(t_observations_used, FBLmean,  'parent',f_axes,'color','r','linewidth',2);

if strcmpi(opts.model, '3s')
    
    statenames = {'c' 'x' 's'};
    
elseif strcmpi(opts.model, '5s2b')
    
    statenames = {'c' 's0' 's1' 's2' 's3' 's4'};
    if ~isfield(moments, 's0_mean')
        
        moments.s0_mean = params.S - moments.s1_mean - moments.s2_mean - moments.s3_mean - moments.s4_mean;
        moments.s0_sd = nan(size(moments.s0_mean));
        
    end
    for k = 1:numel(params.Btot)
        
        statenames{end + 1} = sprintf('b%d', k - 1);
        
    end
    
end
if numel(state_axes) == 1
    
    state_axes = repmat(state_axes, 1, numel(statenames));
    
end
co = get(state_axes(1), 'colororder');
for k = 1:numel(statenames)
    
    j = mod(k - 1, size(co, 1)) + 1;
    y = moments.(sprintf('%s_mean', statenames{k}));
    y_sd = moments.(sprintf('%s_sd', statenames{k}));
    y = y(:); y_sd = y_sd(:);
    
    if normalize_states
        
        y_sd = y_sd / (max(y) - min(y));
        y = (y - min(y)) / (max(y) - min(y));        
        
    end    
    
    if any(~isnan(y_sd))
       
        patchx = reshape([t_states(:); flipud(t_states(:)); t_states(1)],[],1);
        patchcolor = 1 - (1 - co(j, :)) / 2;  % lighten
        patch(patchx, [y - y_sd; flipud(y + y_sd); y(1) - y_sd(1)], patchcolor, 'parent', state_axes(k), 'edgecolor', 'none');
        
    end    
        
    line(t_states, y, 'color', co(j, :), 'parent', state_axes(k));    
    
end
if normalize_states
    
    set(state_axes, 'ylim', [-0.05 1.05]);
    
else

    set(state_axes, 'ylimmode', 'auto');

end

ch = get(f_axes,'children');

set(f_axes,'children',[Fmeanline; fblline; fblpatch; setdiff(ch,[Fmeanline Fintpatch fblline fblpatch]'); Fintpatch]);