function [U_resamp, randperm_resample] = resampling_random_vars(V, T, randtype, gpumode)
%random variables used for resampling
if ~exist('randtype', 'var')
    randtype = 'double';
end
if ~exist('gpumode','var')
    gpumode = false;
end
switch V.resamplemode
    case 'stratified'
        if gpumode
            resamp_variates = gpuArray.rand(V.Nparticles, T, randtype);
        else
            resamp_variates = rand(V.Nparticles, T, randtype);
        end
    case 'systematic'
        if gpumode
            resamp_variates = gpuArray.rand(1, T, randtype);
        else            
            resamp_variates = rand(1, T,randtype);
        end
    otherwise
        error('unrecognized resampling mode');
end
if gpumode
    basevals = gpuArray.colon(cast(0,randtype), cast(V.Nparticles - 1, randtype))';
else
    basevals = cast(0:V.Nparticles - 1, randtype)';
end
if T == 1
    U_resamp = (basevals + resamp_variates) / V.Nparticles;     
else
    U_resamp = bsxfun(@plus, basevals, resamp_variates) / V.Nparticles; 
end
if V.randperm_resample    
    assert(V.Nparticles < 2 ^ 32 - 1, 'too many particles for random permutation');
    randperm_resample = cast((1:V.Nparticles)'*ones(1,T), 'uint32');
    randperm_resample = randperm_jansimon(randperm_resample);
    if gpumode
        randperm_resample = gpuArray(randperm_resampled); %not optimal, ideally we should run the algorithm on the GPU FIXME
    end
else
    randperm_resample = [];
end