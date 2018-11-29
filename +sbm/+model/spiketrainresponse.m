function [brightness, states, tt_states, tt_brightness] = spiketrainresponse(params,opts,dt,apcounts)
T = numel(apcounts);
nsteps = ceil(dt / opts.substepsize); %10 ms increments or so
stepsize = dt / nsteps;

if strcmpi(opts.model, '5s2b')
    
    koff_B = 1 ./ params.tau_B;
    kon_B  = koff_B ./ params.kd_B;
    if ~isfield(params, 'kd_ex'), params.kd_ex = []; end
    
    [seq, beq] = sbm.model.equilibriumstates(params, opts);
    c = nan(1,              T);
    s = nan(numel(seq) + 1, T);
    b = nan(numel(beq),     T);
    c(1) = params.c0;
    s(:, 1) = [params.S - sum(seq); seq];
    b(:, 1) = beq;
    FBGp1 = 1 + params.R * (1 + params.dbrightness * seq / params.S);
    
    for t = 1:T
        
        prevind = max(1, t - 1);
        
        [c(t), nexts, nextb] = sbm.model.simulate(opts.n_newton_iterations, c(prevind) + params.A * apcounts(t), s(:, prevind)', b(:, prevind)', stepsize, params.c0, 1 / params.tau_ex, params.kd_ex, params.kon, params.koff, kon_B, koff_B, params.Btot);
        s(:, t) = nexts';
        b(:, t) = nextb';
        
    end
    
    brightness = params.dbrightness(:)' * s(2:end, nsteps:nsteps:end) / params.S + FBGp1;
    states = [c; s; b];
    
else
    
    error('unrecognized model: %s', opts.model);
    
end

tt_states = (1:T) * stepsize;
tt_brightness = tt_states(nsteps:nsteps:end);