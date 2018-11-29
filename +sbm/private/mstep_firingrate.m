function params = mstep_firingrate(moments, params, opts, dt)
if isempty(params), return; end

dt = reshape(dt, [], 1);
nsteps = ceil(dt / opts.substepsize);
stepsize = reshape(dt ./ nsteps, 1, []); %row vector

spikes_per_segment   = arrayfun(@(v) sum(v.n_mean),   moments); %posterior (smoothed) expectation
substeps_per_segment = arrayfun(@(v) numel(v.n_mean), moments);

if any(isnan(params(1).gamma_kt_FR)) %no prior
    
    fr = sum(spikes_per_segment) / sum(substeps_per_segment .* stepsize); %assigning the posterior mean is the M-step of the EM algorithm for the firing rate
    
else    
    
    [dtlist, ~, dtind] = unique(dt);
    nsteps_unique = ceil(dtlist / opts.substepsize);
    stepsize_unique = reshape(dtlist ./ nsteps_unique, [], 1); %column vector
    
    [W, Wc] = deal(nan(1, numel(dtlist))); %row vectors
    for k = 1:numel(dtlist)
        
        ii = dtind == k;
        W(k)  = sum(  spikes_per_segment(ii));        %sum of weights for all spiking     particles, for this dt
        Wc(k) = sum(substeps_per_segment(ii)) - W(k); %sum of weights for all non-spiking particles, for this dt
        
    end
    
    Qf = @(r) -mstep_firingrate_eLL(r, W, Wc, stepsize_unique, params(1).gamma_kt_FR(1), params(1).gamma_kt_FR(2));
    maxfr = min(1 / max(stepsize));
    
    [fr,fv,ef,op] = fminbnd(Qf, 0, maxfr);

end

for u = 1:numel(params)
    params(u).fr = fr;
end


function Q = mstep_firingrate_eLL(r, W, Wc, stepsize, k, theta)
%ignoring additive offset of -(gammaln(k) + k * log(theta))
Q = W * log(r * stepsize) + Wc * log(1 - r * stepsize) + (k - 1) * log(r) - r / theta;