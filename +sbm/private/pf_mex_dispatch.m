function [marglik, log_sum_raw_w, neff, nmean, gpmean, gpsqmean] = pf_mex_dispatch(V, g, f, u, q_spike0, seedval)

if strcmpi(V.model, '5s2b')
        
    pffunc = @cuda5s_mex;
    
elseif strcmpi(V.model, '3s')
        
    pffunc = @cuda3s_mex;
    
else
    
    error('unrecognized model: %s', V.model);
    
end

if exist('seedval', 'var')
    
    pffunc2 = @(g, f, u, q) pffunc(g, f, u, q, seedval);
    
else
    
    pffunc2 = pffunc;
    
end

if nargout > 3

    [marglik, log_sum_raw_w, neff, nmean, gpmean, gpsqmean] = pffunc2(g, single(f), u, q_spike0);
    
else  % don't compute moments
    
    [marglik, log_sum_raw_w, neff] = pffunc2(g, single(f), u, q_spike0);
    
end