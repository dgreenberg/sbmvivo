function [lik, moments] = pf_mex_wrapper(f,dt,V,P,nA2D,g,u)
min_prewin = 2.5;
if ~exist('nA2D','var'), nA2D = 1; end
if ~exist('u','var'), u = rand(size(f), 'single'); end
T = numel(f);
lik = empty_likstruct;
K = ceil(V.filtersmoother_window / dt);
nsteps = ceil(dt / V.substepsize);
stepsize = dt / nsteps;
ntimepoints_pre = ceil(max(V.filtersmoother_window, min_prewin) / dt);
T_total = T + ntimepoints_pre;
p_spike = min(0.5, P.fr * dt / nsteps);
min_q_spike = 0.01; %1 - exp(log(1 - 0.01) / nsteps);
min_q_spike2 = 0;
q_spike0 = single(ones(1, T_total * nsteps) * max(p_spike, min_q_spike));
gksd = 0.1 / stepsize; %100 ms s.d.
gk_npts = ceil(gksd * 4);
gk = exp(-0.5 * (-gk_npts:gk_npts).^2 / gksd ^ 2); gk = gk / sum(gk); %should deal with edge effects FIXME

if ~exist('g','var')
    
    ginit = true;
    g = pf_mex_init_wrapper(numel(f), V, dt); %also randomly initializes rng seed
    
else
    
    ginit = false; %FIXME check that g is valid?
    
end

if ~V.multiroundpf
    
    g = pf_mex_configure(g, V, P, dt, nA2D, T, nsteps, ntimepoints_pre, K);
    
    if nargout > 1
        
        [marglik, log_sum_raw_w, neff, nmean, gpmean, gpsqmean] = pf_mex_dispatch(V, g, single(f), u, q_spike0);
        
    else
        
        [marglik, log_sum_raw_w, neff] = pf_mex_dispatch(V, g, single(f), u, q_spike0); %don't calculate moments
        
    end
    
else
    
    g = pf_mex_configure(g, V, P, dt, nA2D, T, nsteps, ntimepoints_pre, K);
    
    [~, log_sum_raw_w1, neff1, nmean1] = pf_mex_dispatch(V, g, f, rand(size(f), 'single'), q_spike0);
    assert(~any(isnan(log_sum_raw_w1)) && ~any(isnan(neff1)), 'invalid result');
    q_spike_gpu = max(conv(nmean1, gk, 'same'), min_q_spike2);
    
    g = pf_mex_configure(g, V, P, dt, nA2D, T, nsteps, ntimepoints_pre, K);
    
    if nargout > 1
        
        [marglik, log_sum_raw_w, neff, nmean, gpmean, gpsqmean] = pf_mex_dispatch(V, g, single(f), u, q_spike_gpu);
        
    else
        
        [marglik, log_sum_raw_w, neff] = pf_mex_dispatch(V, g, single(f), u, q_spike_gpu); %don't calculate moments
        
    end
    
end
assert(~any(isnan(log_sum_raw_w)) && all(isfinite(log_sum_raw_w)) && ~any(isnan(neff)), 'invalid result');

lik.logmarglik = marglik;
lik.logcondpnextobs = double(log_sum_raw_w);
lik.neff_proposal = neff;
if ginit
    
    pf_mex_clear(g, V); %clean up
    
end
if nargout > 1
    
    moments = empty_momentstruct(numel(f), nsteps, V);
    nmean = max(0, min(1, double(nmean(nsteps * ntimepoints_pre + 1:end))));
    moments.n_mean = nmean;
    moments.n_sd = sqrt(nmean .* (1 - nmean));
    moments.gp_mean = gpmean;
    moments.gp_sd   = sqrt(gpsqmean - gpmean .^2);
    
end