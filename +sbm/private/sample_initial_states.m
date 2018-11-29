function [c, s, b, n, r, gp, expr, log_w, log_w_corrected, vFtotal, FBGp1] = sample_initial_states(...
    params, opts, F, dt, true_st, nA2D, q_spike_pre)
[Fadj, vF, ~, randtype, spikesknown, ~, stepsize, true_n_sub, p_spike, nrand_r, ~, randseed, ~, ~,~, ~, ...
    ~, ~, ~, ~, ~] ...
    = pfilter_precalc(params, opts, dt, F, true_st, nA2D);

[seq, beq] = sbm.model.equilibriumstates(params, opts);  % calculate equilibrium states
FBGp1 = 1 + params.R * (1 + params.dbrightness * seq / params.S);
seq = [params.S - sum(seq); seq];

nsteps_pre = numel(q_spike_pre);

if opts.singleprecision
    
    urand = nan(opts.Nparticles, nsteps_pre + 1, randtype);
    randseed = quickrand(urand, randseed);
    
else
    
    urand = rand(opts.Nparticles, nsteps_pre + 1, randtype);
    
end

[c, n] = deal(nan(opts.Nparticles, nsteps_pre, randtype));
s = nan(opts.Nparticles, numel(seq), nsteps_pre);
b = nan(opts.Nparticles, numel(beq), nsteps_pre);

n(:) = 0;
if spikesknown
    
    tt_pre = (-(nsteps_pre + 1):0)' * stepsize;
    true_n_sub = histc(true_st, tt_pre);
    
end

c(:, 1) = params.c0;
s(:, :, 1) = repmat(seq(:)', opts.Nparticles, 1);
b(:, :, 1) = repmat(beq(:)', opts.Nparticles, 1);

log_pq_spiking = zeros(opts.Nparticles, 1, randtype);
log_pq_spike = log(p_spike) - log(q_spike_pre);
log_pq_nospike = log(1 - p_spike) - log(1 - q_spike_pre);

%now run the dynamics randomly to generate potential (c,x,s) state configurations prior to the first fluorescence measurment
%FIXME have one common funtion (preferably mex) for this and the pf main loop, and possibly for simulation/fitting as well
koff_B = 1 ./ params.tau_B;
kon_B  = koff_B ./ params.kd_B;

for ii = 1:nsteps_pre
    iiprev = max(1, ii - 1);
    
    if spikesknown
        n(:, ii) = true_n_sub(ii);
    else        
        n(:, ii) = 0;
        spikeind = urand(:,ii) < q_spike_pre(ii);
        log_pq_spiking(spikeind) = log_pq_spiking(spikeind) + log_pq_spike(ii);
        log_pq_spiking(~spikeind) = log_pq_spiking(~spikeind) + log_pq_nospike(ii);
        n(spikeind, ii) = 1; %could do a different proposal here, or a second pass. FIXME?
    end        
        
    [c(:, ii), s(:, :, ii), b(:, :, ii)] = sbm.model.simulate(opts.n_newton_iterations, c(:, iiprev) + n(:, ii) * params.A, s(:, :, iiprev), b(:, :, iiprev), stepsize, params.c0, 1 / params.tau_ex, params.kd_ex, params.kon, params.koff, kon_B, koff_B, params.Btot);    
    
end

brightness = s(:, 2:end, end) * params.dbrightness(:) / params.S + FBGp1;  % brightness given states for each particle
brightness_eq = reshape(seq(2:end), 1, []) * params.dbrightness(:) / params.S + FBGp1;

[r, expr, log_w, log_w_corrected, gp, vFtotal] = sample_initial_r(params, opts, nA2D, F, Fadj, brightness, brightness_eq, nrand_r, randtype, log_pq_spiking, vF, randseed);

if strcmpi(randtype,'single')
    numvars = {'c' 's' 'b' 'n' 'r' 'gp' 'expr' 'log_w' 'log_w_corrected' 'vFtotal'};
    for k = 1:numel(numvars)
        eval([numvars{k} ' = single(' numvars{k} ');']);
    end
end
assert(~any(isnan(log_w)),'NaN weight');