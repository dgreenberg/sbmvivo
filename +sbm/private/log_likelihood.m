function [L, lik, P, logpriorp, moments] = log_likelihood(f, dt, V, P, st, nA2D, gpu3s)
P = reshape(P, 1, []);
calcmoments = nargout > 4;
L = 0;

if V.profile_firing_rate || V.profile_zeta
    
    %FIXME use try/ catch here
    if calcmoments
        
        [L, lik, FR, P, moments] = profile_likelihood(f, dt, V, P, nA2D, st, gpu3s);
        
    else
        
        [L, lik, FR, P] = profile_likelihood(f, dt, V, P, nA2D, st, gpu3s);
        
    end
    
else
    
    for u = 1:numel(P)
        try
            
            FR = P(1).fr;
            
            if V.usegpu
                
                if calcmoments
                    
                    [lik(u), moments(u)] = pf_mex_wrapper(f{u},dt(u),V,P(u),nA2D(u),gpu3s);
                    
                else
                    
                    lik(u) = pf_mex_wrapper(f{u},dt(u),V,P(u),nA2D(u),gpu3s);
                    
                end
                
            else
                
                Vpre = V;
                Vpre.Nparticles = V.Nparticles_prerun;
                Vpre.ResampleThreshold = ceil(V.ResampleThreshold * V.Nparticles_prerun / V.Nparticles);
                [lik1, moments1, moments_pre1] = run_matlab_particlefilter(f{u},dt(u),Vpre,P(u),st{u},nA2D(u));
                
                nsteps = ceil(dt(u) / V.substepsize); %10 ms increments or so
                stepsize = dt(u) / nsteps;
                gksd = 0.1 / stepsize;
                npts = ceil(gksd * 4);
                gk = exp(-0.5 * (-npts:npts).^2 / gksd^2); gk = gk / sum(gk); %should deal with edge effects FIXME
                q_spike = conv([moments_pre1.n_mean moments1.n_mean], gk, 'same'); %should also determine q_spike for time points before first fluorescence measurement FIXME
                if calcmoments
                    
                    [lik(u), moments(u)] = run_matlab_particlefilter(f{u},dt(u),V,P(u),st{u},nA2D(u),[],q_spike);
                    
                else
                    
                    lik(u) = run_matlab_particlefilter(f{u},dt(u),V,P(u),st{u},nA2D(u),[],q_spike);
                    
                end
                
            end
            
        catch ex
            
            fprintf('particle filter failed: %s', ex.message);
            lik(u) = empty_likstruct;
            lik(u).logmarglik = -inf;
            
        end
        
        L = L + lik(u).logmarglik;
        
    end
    
end

%calculate prior probability for these parameters
logpriorp = 0;
if ~any(isnan(P(1).gamma_kt_FR))
    
    gk = P(1).gamma_kt_FR(1);
    gt = P(1).gamma_kt_FR(2);
    
    log_priorp_FR = (gk - 1) * log(FR) - FR / gt - gammaln(gk) - gk * log(gt);
    
    logpriorp = logpriorp + log_priorp_FR;
    
end

if ~any(isnan(P(1).normal_meancov_logSR(:)))
    
    logSRvec = log([P(1).S P(1).R]);
    logpriorp = logpriorp + log(mvnpdf(logSRvec, P(1).normal_meancov_logSR(:, 1)', P(1).normal_meancov_logSR(:, 2:3)));
    
end
assert(~isnan(logpriorp) && ~isnan(L), 'invalid result');