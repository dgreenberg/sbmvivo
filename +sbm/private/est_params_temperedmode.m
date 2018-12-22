function [P, prevP, prevlik] = est_params_temperedmode(opts, params0, dt, f, nA2D)
st = [];
target_mc_LL_variance = 0.1;
ntrials_initial = 15;

params0 = reshape(params0, 1, []);
nseg = numel(params0);
T = cellfun(@numel, f);
assert(numel(T) == nseg && numel(dt) == nseg, 'sizes of dt, P and and f must be the same');
assert(isfield(params0, 'normal_meancov_logSR') && ~any(isnan(params0(1).normal_meancov_logSR(:))), 'lognormal prior on (S, R) must be available');
for u = 1:numel(params0)
    
    params0(u).S = max(opts.absolute_min_S, params0(u).S);
    assert(~opts.shotnoise || params0(u).gain > 0,'initial gain should not be zero in a shotnoise model');
    assert(~opts.darknoise || params0(u).zeta > 0,'initial dark variance should not be zero in a dark noise model');
    
end
if opts.usegpu
    
    gpudata = pf_mex_init_wrapper(T, opts, dt); %randomly sets seedval for gpu rngstates
    
else
    
    gpudata = [];
    
end

%estimate zeta once before R/S
if opts.profile_zeta || opts.profile_firing_rate
    
    [~, ~, ~, params0] = profile_likelihood(f, dt, opts, params0, nA2D, st);
        
end

%get the parameterization, in particular the function PFunc that maps from a vector of unique parameters pvec to the structure array of all parameters for each segment P
[Pfunc, p0, paramnames, index] = paramest_parameterization(params0, opts);

%measure MC variance
if opts.usefiguregraphics, wb = waitbar(0, 'Measuring MC variance'); end
for j = 1:ntrials_initial
    
    try
        if opts.usefiguregraphics, waitbar(j / ntrials_initial, wb); end
    catch
    end
    [Linit(j), ~, ~, priorpinit(j)] = log_likelihood(f, dt, opts, params0, st, nA2D, gpudata);
    
end
if opts.usefiguregraphics && isvalid(wb); close(wb); end

log_posterior_init = Linit + priorpinit; %missing a normalization factor but doesn't matter
temperature = max(std(log_posterior_init) / sqrt(target_mc_LL_variance), opts.tm_mintemp);
fprintf('\n\ninitial MC std.: %f, tempering distribution with temperature %f\n', std(log_posterior_init), temperature);

targetlogpfunc = @(p) LLwrapper(f, dt, opts, Pfunc(p), st, nA2D, gpudata);
proposal = orderfields(struct( ...
    'm', params0(1).normal_meancov_logSR(:, 1)',  ...
    'V', params0(1).normal_meancov_logSR(:, 2:3), ...
    'p', 1 ...    
    ));
box_constraints = [log(opts.absolute_min_S) -inf; inf inf];

displayfunc = [];
if opts.usefiguregraphics
    
    figh = figure('menubar','none', 'units','normalized', 'position',[.1 .1 .8 .8]);    
    
    displayfunc = @(samples, targetlogpvals, acceptance_ratio, accepted, gstar, next_draw, extradata, preliminary_phase, totalt_sec, maxdraws) ...
        tm_aimh_displayfunc(samples, targetlogpvals, acceptance_ratio, accepted, gstar, next_draw, extradata, preliminary_phase, totalt_sec, maxdraws, figh);
    
end

%run the Giordani and Kohn adaptive independent Metropolis-Hasting algorithm
[samples, targetlogpvals, acceptance_ratio, accepted, gstar, extradata, first_update, preliminary_phase_end] = aimh_gk(targetlogpfunc, proposal, p0, temperature, opts.min_iter_paramest, [], [], displayfunc, box_constraints);

if ~isempty(first_update) && size(samples, 1) - first_update >= 100 && sum(accepted(first_update + 1:end)) > 20
    
    gstar = gmfit_gk(samples(first_update + 1:end, :)); %discard the samples before our first Gaussian mixture fit as a "burn-in" (though the justification is somewhat different here than in RWMH)
    
end

p = gmix_ML(gstar); %ML estimate from a Gaussian mixture fit using K-harmonic means, a special initialization techinique and BIC
[~, extradata_final] = targetlogpfunc(p);
P = extradata_final.P;

extradata = cat(2, extradata{:});
prevP = cat(1, extradata.P); %niter x nseg matrix
prevP = cat(1, params0, prevP);
prevlik = cat(1, extradata.lik); %niter x nseg matrix
prevlik = cat(1, prevlik(1,:), prevlik); %hack, FIXME

if opts.usefiguregraphics && isvalid(figh)
    
    close(figh);
    
end

if opts.usegpu
    
    pf_mex_clear(gpudata, opts);
    
end


function [logposterior, extradata] = LLwrapper(f, dt, V, P0, st, nA2D, gpu3s)
%we use this function for compatibility with aimh_gk.m
[LL, lik, P, logpriorp] = log_likelihood(f, dt, V, P0, st, nA2D, gpu3s); %FIXME use one of the inits
logposterior = LL + logpriorp;
extradata = orderfields(struct('LL',LL,'lik',lik,'P',P,'logpriorp',logpriorp));


function tm_aimh_displayfunc(samples, targetlogpvals, acceptance_ratio, accepted, gstar, next_draw, extradata, preliminary_phase, totalt_sec, maxdraws, figh)

if ~isvalid(figh)
    return
end

extradata = cat(2, extradata{:});
allP = cat(1, extradata.P);
S = cat(1, allP(:,1).S);
R = cat(1, allP(:,1).R);

tt = 2 * pi * (0:100) / 100;
dxy = [cos(tt); sin(tt)];

ax = subplot(5,1,1,'parent',figh);
delete(get(ax,'children'));
plot(S, R, 'b.', 'parent',ax);
set(ax,'nextplot','add');
if ~isempty(gstar)    
    %plot the mixture model we've fit to the samples
    plot(exp(gstar.m(:,1)), exp(gstar.m(:,2)),'c+', 'parent',ax);
    for j = 1:numel(gstar.p)
        text(exp(gstar.m(j, 1)), exp(gstar.m(j, 2)), num2str(gstar.p(j)), 'parent',ax); %mixture component weight
        [v,d] = eig(gstar.V(:,:,j));
        isocline_xy = bsxfun(@plus, gstar.m(j, :)', v * sqrt(d) * dxy);
        plot(exp(isocline_xy(1,:)'), exp(isocline_xy(2,:)') , 'c', 'parent',ax);
    end    
    
end
if ~accepted(end)
    plot(exp(next_draw(1)), exp(next_draw(2)), 'rx', 'parent',ax);    
end
plot(S(end), R(end), 'g.', 'parent',ax);
xlabel('S'); ylabel('R');
set(ax,'xlimmode','auto','ylimmode','auto');
if preliminary_phase
    title(ax, sprintf('%d draws, %d accepted, last AR %f (prelim.)', numel(accepted), sum(accepted), acceptance_ratio(end)));
else
    title(ax, sprintf('%d draws, %d accepted, last AR %f', numel(accepted), sum(accepted), acceptance_ratio(end)));
end

ax = subplot(5,1,2,'parent',figh);
plot(S, 'parent',ax);
ylabel('S');

ax = subplot(5,1,3,'parent',figh);
plot(R, 'parent',ax);
ylabel('R');

ax = subplot(5,1,4,'parent',figh);
plot(targetlogpvals, 'parent',ax);
ylabel('log(P)');

ax = subplot(5,1,5,'parent',figh);
plot(cumsum(accepted(:)) ./ (1:numel(accepted))', 'parent',ax);
ylabel('Total acceptance');

time_remaining = (maxdraws - numel(accepted)) * totalt_sec / numel(accepted);

set(figh,'name',sprintf('elapsed time %f s, remaining %f s', totalt_sec, time_remaining));

drawnow;