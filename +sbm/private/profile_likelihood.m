function [loglik, lik, FR, P, moments] = profile_likelihood(f, dt, V, P, nA2D, st, gpudata)
%[loglik, P, moments, likstruct] = profile_likelihood(f, dt, V, P, true_n, nA2D, st, gpudata)
%profiles over firing rates
if ~isa(f, 'cell'), f = {f}; end
if ~exist('st', 'var') || isempty(st), st = cell(size(f)); end

T = cellfun(@numel, f);
nseg = numel(f);
assert(numel(P) == numel(f) && numel(dt) == numel(f) && numel(nA2D) == numel(f), 'P, nA2D and dt must have one element per data segment');

if any(cat(2, P.fr) ~= P(1).fr)
    
    newfr = mean(cat(2, P.fr));
    
    for u = 1:nseg
        
        P(u).fr = newfr;
        
    end
    
end

assert(all(cat(2, P.fr) == P(1).fr), 'firing rate must be constant');
calcmoments = nargout > 4;

min_q_spike_init = 0.01; %1 - exp(log(1 - 0.01) / nsteps);
min_q_spike = 0;

nsteps = ceil(dt / V.substepsize); %vector
stepsize = dt ./ nsteps;  %vector
ntimepoints_pre = ceil(max(V.filtersmoother_window, V.min_prewin) ./ dt); %vector
K = ceil(V.filtersmoother_window ./ dt); %vector
gksd = 0.1 ./ stepsize; %100 ms s.d. vector
gk_npts = ceil(gksd * 4); %vector
T_total = T + ntimepoints_pre;

if V.usegpu && (~exist('gpudata','var') || isempty(gpudata))
    
    ginit = true;
    gpudata = pf_mex_init_wrapper(T, V, dt); %also randomly initializes rng seed
    
else
    
    ginit = false; %FIXME check that g is valid?
    if ~V.usegpu
        
        gpudata = [];
        
    end
    
end

p_spike = min(0.5, max(V.min_firing_rate, P(1).fr) * dt ./ nsteps);

%initialize proposals and filters
[q, gk] = deal(cell(1, nseg));
for u = 1:nseg
    
    gk{u} = exp(-0.5 * (-gk_npts(u):gk_npts(u)) .^ 2 / gksd(u) ^ 2);
    gk{u} = gk{u} / sum(gk{u}); %should deal with edge effects FIXME
    q{u} = single(ones(1, T_total(u) * nsteps(u)) * max(p_spike(u), min_q_spike_init));
    
end

npi = V.n_profile_iterations;
if ~(V.profile_firing_rate || V.profile_zeta)
    
    %nothing to profile
    if V.multiroundpf
        
        npi = 1;
        
    else
        
        npi = 0;
        
    end
    
end

for kk = 1:npi
    
    %get current likelihoods and moments, using a small number of particles
    [L(kk, :), nmean, lik, moments] = profile_likelihood_iteration(V, P, V.Nparticles_prerun, f, dt, nA2D, T, nsteps, ntimepoints_pre, K, q, gpudata, st);
    assert(~any(isnan(L(kk, :))), 'invalid result');
    
    %update proposal distribution on spiking to match current moments, smoothed with a gaussian filter in time.
    for u = 1:nseg
        
        q{u} = max(conv(nmean{u}, gk{u}, 'same'), min_q_spike); %returns single precision array
        
    end
    
    if V.profile_firing_rate
        
        P = mstep_firingrate(moments, P, V, dt); %update the firing rate, over which we are currently profiling
        
    end
    
    if V.profile_zeta
        
        P = mstep_obs(f, moments, nA2D, P, V);
        
    end
    
end

%now do one more estimation with parameters fixed
if calcmoments
    
    [L(npi + 1, :), nmean, lik, moments] = profile_likelihood_iteration(V, P, V.Nparticles, f, dt, nA2D, T, nsteps, ntimepoints_pre, K, q, gpudata, st);
    
else
    
    [L(npi + 1, :), nmean, lik] = profile_likelihood_iteration(V, P, V.Nparticles, f, dt, nA2D, T, nsteps, ntimepoints_pre, K, q, gpudata, st);
    
end
assert(~any(isnan(L(npi + 1, :))), 'invalid result');
%don't update the firing rate again, as we want the final parameter value to match the final likelihood value

loglik = sum(L(npi + 1, :)); %marginal likelihood

if ginit %clear gpu arrays if we initialized them in this function
    
    pf_mex_clear(gpudata, V);
    
end

FR = P(1).fr;


function [L, nmean, lik, moments] = profile_likelihood_iteration(V, P, Nparticles, f, dt, nA2D, T, nsteps, ntimepoints_pre, K, q_spike, g, st)
calcmoments = nargout > 3;
nseg = numel(f);
L = nan(1, nseg);
lik = repmat(empty_likstruct, 1, nseg);
nmean = cell(1, nseg);

%change from the particle count specified in V to some other particle count
V.ResampleThreshold = ceil(V.ResampleThreshold * Nparticles / V.Nparticles);
V.Nparticles = Nparticles;

for u = 1:nseg
    
    if V.usegpu
                
        g = pf_mex_configure(g, V, P(u), dt(u), nA2D(u), T(u), nsteps(u), ntimepoints_pre(u), K(u));
        
        if calcmoments
            
            [L(u), log_sum_raw_w, neff, nmean{u}, gpmean, gpsqmean] = pf_mex_dispatch(V, g, f{u}, rand(size(f{u}), 'single'), q_spike{u});
            
            moments(u) = empty_momentstruct(numel(f{u}), nsteps(u), V);
            moments(u).n_mean = max(0, min(1, double(nmean{u}(nsteps(u) * ntimepoints_pre(u) + 1:end))));
            moments(u).n_sd = sqrt(moments(u).n_mean .* (1 - moments(u).n_mean));
            moments(u).gp_mean = gpmean;
            moments(u).gp_sd   = sqrt(gpsqmean - gpmean .^2);
            
        else
            
            [L(u), log_sum_raw_w, neff] = pf_mex_dispatch(V, g, f{u}, rand(size(f{u}), 'single'), q_spike{u});
            
        end
        
        lik(u).logmarglik = sum(log_sum_raw_w);
        lik(u).logcondpnextobs = double(log_sum_raw_w);
        lik(u).neff_proposal = neff;
        
    else
        
        [lik(u), moments(u), moments_pre] = run_matlab_particlefilter(f{u}, dt(u), V, P(u), st{u}, nA2D(u), [], q_spike{u});        
        L(u) = lik.logmarglik;
        nmean{u} = [moments_pre.n_mean moments(u).n_mean];
        
    end
    
end