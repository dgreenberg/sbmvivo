function [lik, moments, moments_pre, pswarm] = run_filtersmoother(F,dt,V,P,true_st,nA2D,K,q_spike,record_pswarm)
if ~exist('nA2D','var'), nA2D = 1; end
if ~exist('true_st','var'), true_st = []; end
if ~exist('q_spike','var') || isempty(q_spike)
    q_spike = [];
end
if ~exist('K','var') || isempty(K)
    K = ceil(V.filtersmoother_window / dt);
end
if ~exist('record_pswarm', 'var') || isempty(record_pswarm)
    record_pswarm = false;
end
calcmoments = nargout > 1 || record_pswarm;

nbindingsteps = numel(P.kon);
nbuffers = numel(P.Btot);

koff_B = 1 ./ P.tau_B;
kon_B  = koff_B ./ P.kd_B;

lik = empty_likstruct;
[Fadj, vF, T, randtype, spikesknown, nsteps, stepsize, true_n_sub, p_spike, ...
    nrand_r, quickresample, randseed, guardfac, baseresamp, parents_noresamp, equal_log_w, ...
    q_spike, log_pq_spike, log_pq_nospike, log_pq_spiking, ntimepoints_pre] ...
    = pfilter_precalc(P, V, dt, F, true_st, nA2D, q_spike);

F = cast(F, randtype);
urand = nan(V.Nparticles, nsteps, randtype);

[lik.neff, lik.neff_proposal, lik.n_ancestors, lik.logcondpnextobs, lik.logcondpnextobs_proposal] = deal(nan(1,T));
lik.resampled = false(1, T);

%initialize states
%note that we actually run all the pre-date time points as well as one one time point with data, i.e. up to and including the first fluorescence measurement
parent = zeros(V.Nparticles, K + 1, 'uint32');
[c, n, gp, expr, r, vFtotal]  = deal(nan(V.Nparticles, (K + 1) * nsteps, randtype));
s = nan(V.Nparticles, nbindingsteps + 1, (K + 1) * nsteps);
b = nan(V.Nparticles, nbuffers, (K + 1) * nsteps);
moments = empty_momentstruct(T, nsteps, V);
moments_pre = empty_momentstruct(ntimepoints_pre, nsteps, V);
[c_pre, s_pre, b_pre, n_pre, r(:, 1), gp(:, 1), expr(:, 1), log_w, log_w_corrected, vFtotal(:, 1), Fbgp1] = sample_initial_states(P, V, F, dt, true_st, nA2D, ...
    q_spike(1:nsteps * (ntimepoints_pre + 1)));

%initialize weights
if V.singleprecision && V.doubleprecisionweights
    
    log_w = double(log_w); log_w_corrected = double(log_w_corrected);
    baseresamp = (0:V.Nparticles - 1)' / V.Nparticles;
    
end
[w,           log_w,           lik.neff_proposal(1), lik.logcondpnextobs_proposal(1)] = normalize_weights(log_w);
[w_corrected, log_w_corrected, lik.neff(1),          lik.logcondpnextobs(1)] = normalize_weights(log_w_corrected);

if record_pswarm
    
    pswarm = orderfields(struct( ...
        'c',  nan(V.Nparticles, (T + ntimepoints_pre) * nsteps) ...
        ,'n', nan(V.Nparticles, (T + ntimepoints_pre) * nsteps) ...
        ,'w_smoothed', nan(V.Nparticles, T) ...
        ,'w_proposalasprior', nan(V.Nparticles, T) ...
        ,'w_filtered', nan(V.Nparticles, T) ...
        ,'parent', nan(V.Nparticles, T) ...
        ,'r', nan(V.Nparticles, T) ...
        ,'gp', nan(V.Nparticles, T) ...
        ,'s', nan(V.Nparticles, nbindingsteps + 1, (T + ntimepoints_pre) * nsteps) ...
        ,'b', nan(V.Nparticles, nbuffers, (T + ntimepoints_pre) * nsteps) ...
        ,'nsteps', nsteps, 'ntimepoints_pre', ntimepoints_pre ...
        ));
    
    pswarm.c(:, 1:(ntimepoints_pre + 1) * nsteps) = c_pre;
    pswarm.n(:, 1:(ntimepoints_pre + 1) * nsteps) = n_pre;
    pswarm.r(:, 1) = r(:, 1);
    pswarm.gp(:, 1) = gp(:, 1);
    pswarm.s(:, :, 1:(ntimepoints_pre + 1) * nsteps) = s_pre;
    pswarm.b(:, :, 1:(ntimepoints_pre + 1) * nsteps) = b_pre;
    pswarm.w_filtered(:, 1) = w_corrected;    
    pswarm.w_proposalasprior(:, 1) = w;
    
else
    
    pswarm = [];
    
end

%separate pre-data-states from the first states actualy conditioned on data
c(:, 1:nsteps) = c_pre(:, end - nsteps + 1:end);
s(:, :, 1:nsteps) = s_pre(:, :, end - nsteps + 1:end);
b(:, :, 1:nsteps) = b_pre(:, :, end - nsteps + 1:end);
n(:, 1:nsteps) = n_pre(:, end - nsteps + 1:end); %for calculating moments
c_pre = c_pre(:, 1:ntimepoints_pre * nsteps);
s_pre = s_pre(:, :, 1:ntimepoints_pre * nsteps);
b_pre = b_pre(:, :, 1:ntimepoints_pre * nsteps);
n_pre = n_pre(:, 1:ntimepoints_pre * nsteps);

ancestors_present = false(V.Nparticles, 1);
ci = 1; %circular index. indicates which column of the arrays c,r, etc. will hold the next state values
cisub = nsteps;
for t = 2:T
    
    ciprev = ci;
    assert(cisub  == ci * nsteps);
    ci     = mod(ci, K + 1) + 1; %loop around the K + 1 columns
    
    if V.resample_to_proposal
        
        neff_next = lik.neff_proposal(t - 1);
        q_RS = w;
        
    else
        
        neff_next = lik.neff(t - 1);
        q_RS = w_corrected;
        
    end
    
    if neff_next < V.ResampleThreshold %resample particles from the previous time step
        
        lik.resampled(t) = true;
        
        if quickresample
            
            cw = cumsum(q_RS);
            cw = cw / cw(end);
            U_resamp = baseresamp + (rand(1, class(baseresamp)) / V.Nparticles); %no need to call mex file for a single random value for speed, though that would improve reproducibility starting from the same RNG seeds
            parent(:,ci) = resample_wcsnormed_uincreasing(cw, U_resamp);
            
        else
            
            [U_resamp, randperm_resample] = resampling_random_vars(V, 1, randtype, false);
            if V.randperm_resample
                
                q_RS = q_RS(randperm_resample);
                cw = cumsum(q_RS);
                cw = cw / cw(end);
                ii = resample_wcsnormed_uincreasing(cw, U_resamp); %this could be insufficient if variables sizes start exceeding 2^32 - 2 elements. I don't see that happening soon though -dg 04.2014
                parent(:,ci) = randperm_resample(ii);

            else
                
                cw = cumsum(q_RS);
                cw = cw / cw(end);
                parent(:,ci) = resample_wcsnormed_uincreasing(cw, U_resamp);
                
            end
            
        end
        
        nextc = c(parent(:, ci), cisub);
        nexts = s(parent(:, ci), :, cisub);
        nextb = b(parent(:, ci), :, cisub);
        rprev = r(parent(:, ci), ciprev);
        
        if V.resample_to_proposal
            
            log_w_corrected = log_w_corrected(parent(:, ci)) - log_w(parent(:, ci)); %adjust importance weights based on the fact that the sampling distribution doesn't match the target distribution
            [~, log_w_corrected] = normalize_weights(log_w_corrected);
            log_w(:) = equal_log_w;
            
        else
            
            log_w = log_w(parent(:, ci)) - log_w_corrected(parent(:, ci));
            [~, log_w] = normalize_weights(log_w);
            log_w_corrected(:) = equal_log_w;
            
        end
        
    else %no resampling needed
        
        parent(:,ci) = parents_noresamp;
        nextc   = c(:, cisub);
        nexts   = s(:, :, cisub);
        nextb   = b(:, :, cisub);        
        rprev   = r(:, ciprev);
        
        %w, log_w, w_corrected, and log_w_corrected all correspond to
        %normalized probability distributions for the previous step's
        %filtered estimates
    end
    
    %propagate states from time t to time t + 1
    if V.singleprecision
        
        randseed = quickrand(urand, randseed);
        
    else
        
        urand = rand(V.Nparticles, nsteps, randtype);
        
    end
    
    next_subindlist = (ntimepoints_pre + t - 1) * nsteps + 1:(ntimepoints_pre + t) * nsteps;
    next_cisublist = (ci - 1) * nsteps + 1:ci * nsteps;
    
    if spikesknown %insert known spiking
        
        n(:,next_cisublist) = ones(V.Nparticles,1,randtype) * true_n_sub(next_subindlist);
        log_pq_spiking(:) = 0;
        
    else %sample random spiking
        
        n(:,next_cisublist) = cast(bsxfun(@lt, urand, q_spike(next_subindlist)),randtype);
        log_pq_spiking = n(:,next_cisublist) * (log_pq_spike(next_subindlist) - log_pq_nospike(next_subindlist))' + sum(log_pq_nospike(next_subindlist));
        
    end
    
    for jj = 1:nsteps %given spiking, follow c, x and s over short time steps using Ralston's 2nd order RK method
        
        cisub = (ci - 1) * nsteps + jj;
        nextc = nextc + n(:, cisub) * P.A; %increase free calcium due to spiking

        [nextc, nexts, nextb] = sbm.model.simulate(V.n_newton_iterations, nextc, nexts, nextb, stepsize, P.c0, 1 / P.tau_ex, P.kd_ex, P.kon, P.koff, kon_B, koff_B, P.Btot);        
        
        c(:, cisub) = nextc;        
        s(:, :, cisub) = nexts;
        b(:, :, cisub) = nextb;
        
    end
    
    if V.singleprecision
        
        randseed = quickrandn(nrand_r, randseed);
        
    else
        
        nrand_r = randn(V.Nparticles, 1, randtype); %normal variates for state r
        
    end
    
    r(:, ci) = rprev + nrand_r * (P.sigma_r * sqrt(dt));

    expr(:,ci) = exp(r(:,ci));
    
    brightness = s(:, 2:end, cisub) * P.dbrightness(:) / P.S + Fbgp1;  % brightness given states for each particle
    
    gp(:,ci) = expr(:, ci) .* brightness; %"gain times photon flux"
    fresidual = Fadj(t) - gp(:, ci); %difference between expected and observed fluorescence for each particle
    
    %we don't correct for the discrepancy between P(spiking) and q(spiking)
    %until we calculate filtered states or the log likelihood of the next observation
    %this avoids excessive resampling due to discrepancy between P and q
    if V.shotnoise
        
        vFtotal(:, ci) = vF   + (P.gain / nA2D) * gp(:, ci);
        log_pobs      = -0.5 * (fresidual .^ 2 ./ vFtotal(:, ci) + log(2 * pi * vFtotal(:, ci)));
        
    else
        
        log_pobs      = -0.5 * (fresidual .^ 2  / vF            + log(2 * pi * vF));
        
    end
    %right now, log_w and log_w_corrected still represent normalized
    %distributions from the previous step. we will now change this by
    %adding log_pobs to incorporate the next observation and adding log_pobs.
    %we also add log(p/q) to compensate between the mismatch between the prior and
    %proposal
    if V.singleprecision && V.doubleprecisionweights
        
        log_pobs = double(log_pobs); log_pq_spiking = double(log_pq_spiking);
        
    end
    log_w           = log_w           + log_pobs;
    log_w_corrected = log_w_corrected + log_pobs + log_pq_spiking;
    [w_corrected, log_w_corrected, lik.neff(t),          lik.logcondpnextobs(t)         ] = normalize_weights(log_w_corrected); %weights for P(states | F, prior) etc.
    [w,           log_w,           lik.neff_proposal(t), lik.logcondpnextobs_proposal(t)] = normalize_weights(log_w);           %weights for P(states | F, proposal) etc.
    
    if record_pswarm
        
        jj = (ci - 1) * nsteps + (1:nsteps);
        vv = (ntimepoints_pre + t - 1) * nsteps + (1:nsteps);
        
        pswarm.c(:, vv) = c(:, jj);
        pswarm.n(:, vv) = n(:, jj);
        pswarm.r(:, t) = r(:, ci);
        pswarm.gp(:, t) = gp(:, ci);
        pswarm.s(:, :, vv) = s(:, :, jj);
        pswarm.b(:, :, vv) = b(:, :, jj);
        pswarm.w_proposalasprior(:, t) = w;
        pswarm.w_filtered(:, t) = w_corrected;        
        pswarm.parent(:, t) = parent(:, ci);
        
    end   
    
    assert(~any(isnan(w)) && ~any(isnan(w_corrected)),'NaN weight');
    if calcmoments && t > K %trace particle ancestries backward in time to find the states we'll use to compute the filter-smoother estimates
        
        aci = ci;
        ancestors = parent(:,ci);
        aci = mod(aci - 2, K + 1) + 1; %loop around the K + 1 columns in the opposite direction, going backward in time
        for j = 2:K %FIXME this loop should be recoded in C
            
            ancestors = parent(ancestors, aci);
            aci = mod(aci - 2, K + 1) + 1; %loop around the K + 1 columns in the opposite direction, going backward in time
            
        end
        
        if record_pswarm
            
            % we need to do a reduction on weights
            pswarm.w_smoothed(:, t - K) = 0;
            for au = unique(ancestors)'            
                
                pswarm.w_smoothed(au, t - K) = sum(w_corrected(ancestors == au));
                
            end
            
        end
        
        ancestors_present(:) = false;
        ancestors_present(ancestors) = true;
        lik.n_ancestors(t - K) = sum(ancestors_present);
        %calculate moments for states with one value per observation:
        moments.r_mean(t - K)       =      w_corrected' *  r(ancestors, aci);
        moments.r_sd(t - K)         = sqrt(w_corrected' * (r(ancestors, aci)    - moments.r_mean(t - K)) .^ 2);
        moments.expr_mean(t - K)    =      w_corrected' *  expr(ancestors, aci);
        moments.expr_sd(t - K)      = sqrt(w_corrected' * (expr(ancestors, aci) - moments.expr_mean(t - K)) .^ 2);
        moments.gp_mean(t - K)      =      w_corrected' *  gp(ancestors, aci);
        moments.gp_sd(t - K)        = sqrt(w_corrected' * (gp(ancestors, aci)   - moments.gp_mean(t - K)) .^ 2);
        if V.shotnoise
            
            moments.Fpred_sd(t - K) = sqrt(moments.gp_sd(t - K) ^ 2 + w_corrected' * vFtotal(ancestors, aci));
            
        else
            
            moments.Fpred_sd(t - K) = sqrt(moments.gp_sd(t - K) ^ 2 + vF);
            
        end
        %calculate moments for states with nsteps values per observation:
        subaci             = (1:nsteps) + nsteps * (aci - 1);   %indices into state value arrays, we need multiple columns for this single observation
        moment_index_range = (1:nsteps) + nsteps * (t - K - 1); %indices of the sub-sampling resolution moments we're about to assign
        moments.n_mean(moment_index_range) =             w_corrected' * n(ancestors, subaci);
        moments.n_sd(moment_index_range)   = sqrt(max(0, w_corrected' * n(ancestors, subaci) .^ 2 - moments.n_mean(moment_index_range) .^ 2));
        moments.c_mean(moment_index_range) =             w_corrected' * c(ancestors, subaci);
        moments.c_sd(moment_index_range)   = sqrt(max(0, w_corrected' * c(ancestors, subaci) .^ 2 - moments.c_mean(moment_index_range) .^ 2));
        
        s1 = reshape(s(ancestors, 2, subaci), V.Nparticles, nsteps);
        moments.s1_mean(moment_index_range) = w_corrected' * s1;
        moments.s1_sd(moment_index_range)   = sqrt(max(0, w_corrected' * s1 .^ 2 - moments.s1_mean(moment_index_range) .^ 2));
        
        s2 = reshape(s(ancestors, 3, subaci), V.Nparticles, nsteps);
        moments.s2_mean(moment_index_range) = w_corrected' * s2;
        moments.s2_sd(moment_index_range)   = sqrt(max(0, w_corrected' * s2 .^ 2 - moments.s2_mean(moment_index_range) .^ 2));
        
        s3 = reshape(s(ancestors, 4, subaci), V.Nparticles, nsteps);
        moments.s3_mean(moment_index_range) = w_corrected' * s3;
        moments.s3_sd(moment_index_range)   = sqrt(max(0, w_corrected' * s3 .^ 2 - moments.s3_mean(moment_index_range) .^ 2));
        
        s4 = reshape(s(ancestors, 5, subaci), V.Nparticles, nsteps);
        moments.s4_mean(moment_index_range) = w_corrected' * s4;
        moments.s4_sd(moment_index_range)   = sqrt(max(0, w_corrected' * s4 .^ 2 - moments.s4_mean(moment_index_range) .^ 2));
        
        s0 = P.S - s1 - s2 - s3 - s4;
        moments.s0_mean(moment_index_range) = w_corrected' * s0;
        moments.s0_sd(moment_index_range)   = sqrt(max(0, w_corrected' * s0 .^ 2 - moments.s0_mean(moment_index_range) .^ 2));
        
        b0 = reshape(b(ancestors, 1, subaci), V.Nparticles, nsteps);
        moments.b0_mean(moment_index_range) = w_corrected' * b0;
        moments.b0_sd(moment_index_range)   = sqrt(max(0, w_corrected' * b0 .^ 2 - moments.b0_mean(moment_index_range) .^ 2));
        
        b1 = reshape(b(ancestors, 2, subaci), V.Nparticles, nsteps);
        moments.b1_mean(moment_index_range) = w_corrected' * b1;
        moments.b1_sd(moment_index_range)   = sqrt(max(0, w_corrected' * b1 .^ 2 - moments.b1_mean(moment_index_range) .^ 2));
        
        %moments of states before first fluorescence measurement
        if t == K + 1
            
            moments_pre.n_mean = w_corrected' * n_pre(ancestors, :);
            moments_pre.n_sd   = sqrt(max(0, w_corrected' * n_pre(ancestors, :) .^ 2 - moments_pre.n_mean .^ 2));
            moments_pre.c_mean = w_corrected' * c_pre(ancestors, :);
            moments_pre.c_sd   = sqrt(max(0, w_corrected' * c_pre(ancestors, :) .^ 2 - moments_pre.c_mean .^ 2));            
            
            s1 = reshape(s_pre(ancestors, 2, :), V.Nparticles, []);
            moments_pre.s1_mean = w_corrected' * s1;
            moments_pre.s1_sd   = sqrt(max(0, w_corrected' * s1 .^ 2 - moments_pre.s1_mean .^ 2));
            
            s2 = reshape(s_pre(ancestors, 3, :), V.Nparticles, []);
            moments_pre.s2_mean = w_corrected' * s2;
            moments_pre.s2_sd   = sqrt(max(0, w_corrected' * s2 .^ 2 - moments_pre.s2_mean .^ 2));
            
            s3 = reshape(s_pre(ancestors, 4, :), V.Nparticles, []);
            moments_pre.s3_mean = w_corrected' * s3;
            moments_pre.s3_sd   = sqrt(max(0, w_corrected' * s3 .^ 2 - moments_pre.s3_mean .^ 2));
            
            s4 = reshape(s_pre(ancestors, 5, :), V.Nparticles, []);
            moments_pre.s4_mean = w_corrected' * s4;
            moments_pre.s4_sd   = sqrt(max(0, w_corrected' * s4 .^ 2 - moments_pre.s4_mean .^ 2));
            
            s0 = P.S - s1 - s2 - s3 - s4;
            moments_pre.s0_mean = w_corrected' * s0;
            moments_pre.s0_sd   = sqrt(max(0, w_corrected' * s0 .^ 2 - moments_pre.s0_mean .^ 2));
            
            b0 = reshape(b_pre(ancestors, 1, :), V.Nparticles, []);
            moments_pre.b0_mean = w_corrected' * b0;
            moments_pre.b0_sd   = sqrt(max(0, w_corrected' * b0 .^ 2 - moments_pre.b0_mean .^ 2));
            
            b1 = reshape(b_pre(ancestors, 2, :), V.Nparticles, []);
            moments_pre.b1_mean = w_corrected' * b1;
            moments_pre.b1_sd   = sqrt(max(0, w_corrected' * b1 .^ 2 - moments_pre.b1_mean .^ 2));
            
        end
        
    end
    
end

lik.logmarglik = sum(lik.logcondpnextobs);

if ~calcmoments, return; end
%compute the smoother states for final K observations where looking forward in time by a full K time steps is not possible so we look to w at T instead
for t = max(1, T - K + 1):T
    
    %note that the current value of w is for t == T
    aci = ciprev; %corresponds to time step T
    if t == T %for the final time step, the filter-smoother is simply the filter
        
        ancestors = uint32((1:V.Nparticles)');
        
    else
        
        ancestors = parent(:,ci);
        aci = mod(aci - 2, K + 1) + 1; %loop around the K + 1 columns in the opposite direction, going backward in time
        for j = 2:T - t
            
            ancestors = parent(ancestors, aci);
            aci = mod(aci - 2, K + 1) + 1; %loop around the K + 1 columns in the opposite direction, going backward in time
            
        end
        
    end
    if record_pswarm
        
        % we need to do a reduction on weights
        pswarm.w_smoothed(:, t) = 0;
        for au = unique(ancestors)'
            
            pswarm.w_smoothed(au, t) = sum(w_corrected(ancestors == au));
            
        end
        
    end
    
    ancestors_present(:) = false;
    ancestors_present(ancestors) = true;
    lik.n_ancestors(t) = sum(ancestors_present);
    %could speed things up a bit here by putting the moment calculation inside the loop over j. FIXME?
    moments.r_mean(t)       =      w_corrected' *  r(ancestors, aci);
    moments.r_sd(t)         = sqrt(w_corrected' * (r(ancestors, aci)    - moments.r_mean(t)) .^ 2);
    moments.expr_mean(t)    =      w_corrected' *  expr(ancestors, aci);
    moments.expr_sd(t)      = sqrt(w_corrected' * (expr(ancestors, aci) - moments.expr_mean(t)) .^ 2);
    moments.gp_mean(t)      =      w_corrected' *  gp(ancestors, aci);
    moments.gp_sd(t)        = sqrt(w_corrected' * (gp(ancestors, aci)   - moments.gp_mean(t)) .^ 2);
    if V.shotnoise
        
        moments.Fpred_sd(t) = sqrt(moments.gp_sd(t) ^ 2 + w_corrected' * vFtotal(ancestors, aci));
        
    else
        
        moments.Fpred_sd(t) = sqrt(moments.gp_sd(t) ^ 2 + vF);
        
    end
    
    %calculate moments for states with nsteps values per observation:
    subaci             = (1:nsteps) + nsteps * (aci - 1);   %indices into state value arrays, we need multiple columns for this single observation
    moment_index_range = (1:nsteps) + nsteps * (t - 1);     %indices of the sub-sampling resolution moments we're about to assign
    moments.n_mean(moment_index_range) =             w_corrected' * n(ancestors, subaci);
    moments.n_sd(moment_index_range)   = sqrt(max(0, w_corrected' * n(ancestors, subaci) .^ 2 - moments.n_mean(moment_index_range) .^ 2));
    moments.c_mean(moment_index_range) =             w_corrected' * c(ancestors, subaci);
    moments.c_sd(moment_index_range)   = sqrt(max(0, w_corrected' * c(ancestors, subaci) .^ 2 - moments.c_mean(moment_index_range) .^ 2));
    
    s1 = reshape(s(ancestors, 2, subaci), V.Nparticles, nsteps);
    moments.s1_mean(moment_index_range) = w_corrected' * s1;
    moments.s1_sd(moment_index_range)   = sqrt(max(0, w_corrected' * s1 .^ 2 - moments.s1_mean(moment_index_range) .^ 2));
    
    s2 = reshape(s(ancestors, 3, subaci), V.Nparticles, nsteps);
    moments.s2_mean(moment_index_range) = w_corrected' * s2;
    moments.s2_sd(moment_index_range)   = sqrt(max(0, w_corrected' * s2 .^ 2 - moments.s2_mean(moment_index_range) .^ 2));
    
    s3 = reshape(s(ancestors, 4, subaci), V.Nparticles, nsteps);
    moments.s3_mean(moment_index_range) = w_corrected' * s3;
    moments.s3_sd(moment_index_range)   = sqrt(max(0, w_corrected' * s3 .^ 2 - moments.s3_mean(moment_index_range) .^ 2));
    
    s4 = reshape(s(ancestors, 5, subaci), V.Nparticles, nsteps);
    moments.s4_mean(moment_index_range) = w_corrected' * s4;
    moments.s4_sd(moment_index_range)   = sqrt(max(0, w_corrected' * s4 .^ 2 - moments.s4_mean(moment_index_range) .^ 2));
    
    s0 = P.S - s1 - s2 - s3 - s4;
    moments.s0_mean(moment_index_range) = w_corrected' * s0;
    moments.s0_sd(moment_index_range)   = sqrt(max(0, w_corrected' * s0 .^ 2 - moments.s0_mean(moment_index_range) .^ 2));
    
    b0 = reshape(b(ancestors, 1, subaci), V.Nparticles, nsteps);
    moments.b0_mean(moment_index_range) = w_corrected' * b0;
    moments.b0_sd(moment_index_range)   = sqrt(max(0, w_corrected' * b0 .^ 2 - moments.b0_mean(moment_index_range) .^ 2));
    
    b1 = reshape(b(ancestors, 2, subaci), V.Nparticles, nsteps);
    moments.b1_mean(moment_index_range) = w_corrected' * b1;
    moments.b1_sd(moment_index_range)   = sqrt(max(0, w_corrected' * b1 .^ 2 - moments.b1_mean(moment_index_range) .^ 2));
        
end
moments.Fpred_mean = moments.gp_mean + P.fdc;