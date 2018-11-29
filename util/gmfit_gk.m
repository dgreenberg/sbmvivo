function gmix = gmfit_gk(data, maxclusters, maxK, exponent, subsample_frac, maxabsskewness)
%gmix = gmfit_gk(data, maxclusters, maxK, subsample_frac, maxabsskewness)
%
%data is nsamples x datadim
%gmix.m is K x datadim array of mixture component means
%gmix.V is datadim x datadim x K array of mixture component covariances
%gmix.p is a 1 x K vector of component weight
%
%Robust and fast fitting of Gaussian mixtures using the procedure describe in
%"Adaptive Independent Metropolis-Hastings by Fast Estimation of Mixtures of Normals" Giordani & Kohn, 2008
%http://amstat.tandfonline.com/doi/abs/10.1198/jcgs.2009.07174
[nsamples, datadim] = size(data);
nsamples_unique = size(unique(data, 'rows'), 1);
assert(nsamples_unique >= 5, 'at least five unique samples required');

if ~exist('exponent', 'var') || isempty(exponent)
    exponent = []; %use default value inside kmeanhar
end
if ~exist('subsample_frac', 'var') || isempty(subsample_frac)
    subsample_frac = 0.1;
end
if ~exist('maxK','var') || isempty(maxK)
    maxclusters = min(ceil(nsamples_unique / 2), 5);
end
if ~exist('maxskewness','var') || isempty(maxskewness)
    maxabsskewness = 0.2;
end

normvalvars = abs(skewness(data)) < maxabsskewness; %FIXME not yet used
covdata = cov(data);

[m_allK, V_allK, p_allk] = deal(cell(1, maxclusters));
LL = nan(1, maxclusters);
for nclusters = 1:maxclusters
    
    cinit = clusterinit(data, nclusters, [], [], subsample_frac);
    [m_allK{nclusters}, ~, ~, ~, qik] = kmeanhar(data, nclusters, [], exponent, cinit);
    assert(~any(isnan(m_allK{nclusters}(:))) && ~any(isnan(qik(:))), 'invalid clustering output');
    %estimate variance and cluster weight for each cluster. note that this is based on the KHM error function and is NOT a Gaussian mixture ML estimate
    %the idea is to make it more robust especially in the presence of rejected draws etc.
    %qik(k,i) is equal to q_ik in [Zhang2000]
    
    qk = sum(qik, 2);
    pik = bsxfun(@rdivide, qik, qk);
    V_allK{nclusters} = nan(datadim, datadim, nclusters);
    for k = 1:nclusters
        
        d = bsxfun(@minus, data, m_allK{nclusters}(k, :)); %displacement of all data from the i-th cluster center
        
        V_allK{nclusters}(:, :, k) = bsxfun(@times, d', pik(k, :)) * d;
        
        if all(all(V_allK{nclusters}(:, :, k) == 0))
            
            V_allK{nclusters}(:, :, k) = eye(datadim); %to avoid nans when calculating likelihood
            
        elseif any(eig(V_allK{nclusters}(:, :, k)) < eps) %ensure covarance matrix is positive definite
            
            V_allK{nclusters}(:, :, k) = 0.25 * covdata;
            
        end
        
    end
    
    p_allk{nclusters} = qk' / sum(qk);
    
    assert(~any(isnan(V_allK{nclusters}(:))) && ~any(isnan(p_allk{nclusters}(:))), 'invalid Gaussian mixture');
    
    %calculate log-likelihood of the data given this model
    LL_all = zeros(nclusters, nsamples);
    for k = 1:nclusters
        
        d = bsxfun(@minus, data, m_allK{nclusters}(k, :))'; %displacement of all data from the i-th cluster center
        dtSinvd = sum((V_allK{nclusters}(:, :, k) \ d) .* d, 1);
        LL_all(k, :) = -0.5 * (dtSinvd + log(det(V_allK{nclusters}(:, :, k)) * (2 * pi) ^ datadim)) + log(p_allk{nclusters}(k));
        
    end
    
    LL_best = max(LL_all, [], 1); %LL from best cluster center for each data point
    
    LL_eachsample = log(sum(exp(bsxfun(@minus, LL_all, LL_best)), 1)) + LL_best; % subtract off max before using exp to avoid numerical issues
    
    LL(nclusters) = sum(LL_eachsample);
    
end

nparameters = (1:maxclusters) * (2 + 3 + 1); %means, variances, covariance, mixture component probability.
BICvals = -2 * LL + nparameters * log(nsamples);
[~, bestnclusters] = min(BICvals);
m = m_allK{bestnclusters};
V = V_allK{bestnclusters};
p = p_allk{bestnclusters};
gmix = orderfields(struct('m', m, 'V', V, 'p', p));