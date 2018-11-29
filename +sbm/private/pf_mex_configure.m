function gpudata = pf_mex_configure(gpudata, V, P, dt, nA2D, T, nsteps, ntimepoints_pre, K)

if strcmpi(V.model, '5s2b')
    
    [seq, ~] = sbm.model.equilibriumstates(P, V);  % calculate equilibrium states
    FBGp1 = 1 + P.R * (1 + P.dbrightness * seq / P.S);
    gpudata = cuda5s_mex_configure(gpudata, V.Nparticles, P, dt, nA2D, T, nsteps, ntimepoints_pre, K, V.ResampleThreshold, V.n_newton_iterations, FBGp1);
    
elseif strcmpi(V.model, '3s')
    
    gpudata = cuda3s_mex_configure(gpudata, V.Nparticles, P, dt, nA2D, T, nsteps, ntimepoints_pre, K, V.ResampleThreshold);
    
else
    
    error('unrecognized model: %s', V.model);
    
end