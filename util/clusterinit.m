function c = clusterinit(data, K, clusterfunc, J, subsample_frac, c0)
%data is nsamples x datadim
%c is K x datadim
%clusterfunc operates as
%[centers, objfunval] = clusterfunc(data, K, centers0)
%
%"Refining Initial Points for K-Means Clustering" Bradley & Fayyad
%http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=2B5A4023CB028CC5400710AAB75EBC84?doi=10.1.1.44.5872&rep=rep1&type=pdf
[nsamples, datadim] = size(data);
assert(nsamples >= K,'need at least as many samples as clusters');

if ~exist('K','var') || isempty(K)
    K = 3;
end
if ~exist('clusterfunc','var') || isempty(clusterfunc)
    clusterfunc = @(data, K, c0) kmeanhar(data, K, [], 2, c0);
end
if ~exist('J','var') || isempty(J)
    J = 10;
end
if ~exist('subsample_frac','var') || isempty(subsample_frac)
    subsample_frac = 0.1;
end
if ~exist('c0','var') || isempty(c0)
    c0 = data(randsubset(nsamples, K), :); %random initialization
end

nsubsamples = ceil(nsamples * subsample_frac);

if nsubsamples < K
    %warning('insufficient samples to initialize clustering based on sub-sampling, using random initialization');
    c = c0;
    return;
end

CM = nan(K * J, datadim);
for i = 1:J
    
    r = randsubset(nsamples, nsubsamples);
    row_range = (1:K) + (i - 1) * K;
    CM(row_range, :) = clusterfunc(data(r, :), K, c0);
    
end

CM_ok = CM(~any(isnan(CM), 2), :);

FM = nan(K * J, datadim);
ov = nan(1, J);
for i = 1:J
    
    row_range = (1:K) + (i - 1) * K;
    [FM(row_range, :), ov(i)] = clusterfunc(CM_ok, K, CM(row_range, :));
    
end

[~, besti] = min(ov);
row_range = (1:K) + (besti - 1) * K;
c = FM(row_range, :);
assert(~any(isnan(c(:))), 'invalid clustering result');


function s = randsubset(nsamples, nsubsamples)
r = randperm(nsamples); %FIXME more efficiency for large nsamples?
s = r(1:nsubsamples);

% function c = clustermod(clusterfunc, data, K, c0)
% c = clusterfunc(data, K, c0);
% bsxfun(@minus, c', data);