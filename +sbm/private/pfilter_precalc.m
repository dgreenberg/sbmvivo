function [Fadj, vF, T, randtype, spikesknown, nsteps, stepsize, true_n_sub, p_spike, ...
    nrand_r, quickresample, randseed, guardfac, baseresamp, parents_noresamp, equal_log_w, ...
    q_spike, log_pq_spike, log_pq_nospike, log_pq_spiking, ntimepoints_pre] ...
    = pfilter_precalc(P, V, dt, F, true_st, nA2D, q_spike)
min_prewin = 2.5;
nsteps = ceil(dt / V.substepsize); %10 ms increments or so
ntimepoints_pre = ceil(max(V.filtersmoother_window, min_prewin) / dt);
stepsize = dt / nsteps;
min_q_spike = 0.01; %1 - exp(log(1 - 0.01) / nsteps);

Fadj = F - P.fdc;
if V.darknoise
    
    vF = P.zeta / nA2D;
    
else
    
    vF = 0;
    
end
T = numel(F);
T_total = T + ntimepoints_pre;
equal_log_w = log(1 / V.Nparticles);
spikesknown = ~isempty(true_st) && V.use_true_n_estep;
quickresample = ~V.randperm_resample && strcmpi(V.resamplemode,'systematic');
if spikesknown
    
    it = (-ntimepoints_pre:numel(F))' * dt;
    tt = reshape(bsxfun(@plus, reshape(it,1,[]) - dt, dt * (1:nsteps)' / nsteps),[],1);    
    true_n_sub = histc(true_st, tt);
    true_n_sub = true_n_sub(:)';
    p_spike = NaN;
    
else
    
    true_n_sub = [];
    p_spike = P.fr * dt / nsteps;
    if p_spike > 0.5
        
        warning('sbm:FRtoohigh','firing rate lowered from %f to maximum of %f Hz for a step size of %f', P.fr, 0.5 / (dt / nsteps), dt / nsteps);
        p_spike = 0.5;
        
    end
    
end
if ~exist('q_spike','var') || isempty(q_spike)
    
    q_spike = ones(1, T_total * nsteps) * max(p_spike, min_q_spike);    
    
elseif numel(q_spike) == 1
    
    q_spike = ones(1, T_total * nsteps) * q_spike;
    
else
    
    assert(numel(q_spike) == T_total * nsteps, 'q_spike is incorrectly sized');
    q_spike = max(q_spike, p_spike);
    
end
log_pq_spike = log(p_spike) - log(q_spike);
log_pq_nospike = log(1 - p_spike) - log(1 - q_spike);
log_pq_spiking = zeros(V.Nparticles, 1); %we will use this to calculate the ratio of p(spiking) / q(spiking) in order to calculate importance weights
nrand_r = nan(V.Nparticles, 1);
randseed = randi(2^32 - 1, 1,'uint32');
guardfac = 1 / 0.98;
if V.singleprecision
    
    randtype = 'single';
    numvars = {'Fadj' 'vF' 'T' 'p_spike' 'stepsize' 'nrand_r' ...
        'guardfac' 'equal_log_w' 'q_spike' 'log_pq_spike' 'log_pq_nospike' 'log_pq_spiking'};
    for k = 1:numel(numvars)
        
        eval([numvars{k} ' = single(' numvars{k} ');']);
        
    end
    
else
    
    randtype = 'double';
    
end
if quickresample
    
    baseresamp = cast((0:V.Nparticles - 1) / V.Nparticles, randtype)';    
    
else
    
    baseresamp = [];
    
end
parents_noresamp = uint32((1:V.Nparticles)');