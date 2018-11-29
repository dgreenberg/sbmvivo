function pf_mex_clear(gpudata, V)

if strcmpi(V.model, '5s2b')
    
    cuda5s_mex_clear(gpudata);
    
elseif strcmpi(V.model, '3s')
    
    cuda3s_mex_clear(gpudata);
    
else
    
    error('unrecognized model: %s', V.model);
    
end