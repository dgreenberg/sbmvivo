function [x, maxlikelihood] = gmix_ML(gmix)

[n, d] = size(gmix.m);

Cinv = nan(d, d, n);
Cdet = nan(1, n);
mCinv = nan(n, d);
for k = 1:n
    
    Cinv(:, :, k) = inv(gmix.V(:, :, k));
    Cdet(1, k)    = det(gmix.V(:, :, k));
    mCinv(k, :)   = gmix.m(k, :) * Cinv(:, :, k);
    
end
alpha = gmix.p(:) ./ sqrt((2 * pi) ^ d * Cdet(:));

ofunc = @(x) -gmix_ML_p_and_gradient(x, gmix, n, d, Cinv, Cdet, mCinv, alpha);
allx = nan(d, n);
allp = nan(1, n);
for k = 1:n
    
    x0 = gmix.m(k, :)';
    [allx(:, k), negp, ef, op] = fminsearch(ofunc, x0, optimset('tolfun', 1e-8, 'tolx', 1e-8)); %FIXME switch to a gradient-based method
    allp(k) = -negp;
    
end

[maxlikelihood, bestk] = max(allp);
x = allx(:, bestk);

function [p, g] = gmix_ML_p_and_gradient(x, gmix, n, d, Cinv, Cdet, mCinv, alpha)
dxt = bsxfun(@minus, x', gmix.m);

p_eachcomp = nan(n, 1);
g_eachcomp = nan(n, d);
for k = 1:n
   
   dxtCinv = dxt(k, :) * Cinv(:, :, k);
   p_eachcomp(k) = exp(-0.5 * dxtCinv * dxt(k, :)') * alpha(k);

   if nargout > 1
   
       g_eachcomp(k, :) = p_eachcomp(k) * (mCinv(k, :) - dxtCinv);
       
   end
   
end

p = sum(p_eachcomp);

if nargout > 1
    
    g = sum(g_eachcomp)';
    
end