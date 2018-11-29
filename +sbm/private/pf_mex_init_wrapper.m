function g = pf_mex_init_wrapper(T, V, dt, seedval)
if ~exist('seedval','var')
    seedval = randi_u64;
else
    assert(strcmpi(class(seedval), 'uint64') && numel(seedval) == 1,'seedval must be a uint64 scalar');
end
maxT = max(T);
nsteps = ceil(dt / V.substepsize);
maxnsteps = max(nsteps);
K = ceil(V.filtersmoother_window ./ dt);
maxK = 2 ^ nextpow2(max(K) + 1) - 1; %need K+1 to be a power of 2 so we can do fast modular arithmetic
maxKp1substeps = max((K + 1) .* nsteps);
ntimepoints_pre = ceil(max(V.filtersmoother_window, V.min_prewin) ./ dt);
maxpresubsteps = max(ntimepoints_pre .* nsteps);
maxtotalsubsteps = max((ntimepoints_pre + T) .* nsteps);
if strcmpi(V.model, '5s2b')
    
    g = cuda5s_mex_init(maxT, max(V.Nparticles, V.Nparticles_prerun), maxnsteps, maxK, maxKp1substeps, maxtotalsubsteps, maxpresubsteps, seedval);
    
elseif strcmpi(V.model, '3s')
    
    g = cuda3s_mex_init(maxT, max(V.Nparticles, V.Nparticles_prerun), maxnsteps, maxK, maxKp1substeps, maxtotalsubsteps, maxpresubsteps, seedval);
    
else
    
    error('unrecognized model: %s', V.model);
    
end