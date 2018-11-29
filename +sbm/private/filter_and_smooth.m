function [moments, lik, M] = filter_and_smooth(F, dt, V, P, nA2D, st)
nseg = numel(F);

lik = repmat(empty_likstruct, size(P));
moments = repmat(empty_momentstruct, size(P));
for u = 1:nseg
    
    if V.verbose > 3 && numel(F) > 1
        fprintf('Segment %d/%d ',u,numel(F));
    end
    [moments(u), lik(u), M(u)] = filter_and_smooth_single_segment(V, P(u), dt(u), F{u}, nA2D(u), st{u});
    
end

if nargout > 2 && strcmpi(V.parameterestimation,'em')
    M = reshape(orderfields(M), size(P));
end

function [moments, lik, nextM] = filter_and_smooth_single_segment(V, P, dt, f, nA2D, st)
if V.usegpu
    
    [lik, moments] = pf_mex_wrapper(f,dt,V,P,nA2D); %fixme should keep from reallocating CUDA memory for each segment
    
else
    
    nsteps = ceil(dt / V.substepsize); %10 ms increments or so
    stepsize = dt / nsteps;
    gksd = 0.1 / stepsize;
    npts = ceil(gksd * 4);
    gk = exp(-0.5 * (-npts:npts).^2 / gksd^2); gk = gk / sum(gk); %should deal with edge effects FIXME
    
    if ~V.multiroundpf
        
        [lik, moments] = run_matlab_particlefilter(f, dt, V, P, st, nA2D);
        
    else
        
        V2 = V;
        V2.Nparticles = V.Nparticles_prerun;
        V2.ResampleThreshold = ceil(V.ResampleThreshold * V.Nparticles_prerun / V.Nparticles);
        [lik1, moments1, moments_pre1] = run_matlab_particlefilter(f,dt,V2,P,st,nA2D);
        q_spike = conv([moments_pre1.n_mean moments1.n_mean], gk, 'same');
        [lik, moments] = run_matlab_particlefilter(f,dt,V,P,st,nA2D,[],q_spike);
        
    end
    
end

nextM = struct();