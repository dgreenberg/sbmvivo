function [lik, moments, moments_pre] = run_matlab_particlefilter(F, dt, V, P, true_st, nA2D, K, q_spike)
if ~exist('q_spike','var') || isempty(q_spike)
    
    q_spike = [];
    
end
if ~exist('K','var') || isempty(K)
    
    K = ceil(V.filtersmoother_window / dt);
    
end

if strcmpi(V.model, '5s2b')
    
    [lik, moments, moments_pre] = run_filtersmoother(F,dt,V,P,true_st,nA2D,K,q_spike);
    
else
    
    error('unrecognized model: %s', V.model);
    
end