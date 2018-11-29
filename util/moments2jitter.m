function [st, stsd, gfit] = moments2jitter(ns, maxjitter)
%[st, stsd, gfit] = moments2jitter(ns, maxjitter)
%fit spike train with jitter to posterior spike count moments

wstate = warning('query', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

sdfac = 3; %big enough so that a window this big contains nearly all the area under the curve. normcdf(3) - normcdf(-3) = 0.9973
minjitter = 0.5;

ns = reshape(ns, [], 1);
if ~exist('maxjitter','var')
    maxjitter = 50;
end
maxjitter = max(maxjitter, minjitter + 0.5);

ntest = 100;
sigma_test_list = minjitter + (maxjitter - minjitter) * ...
    (0.01 + 0.98 * ((0:(ntest - 1)) / (ntest - 1)) ); %test values range from 1% to 99% between min and max
T = numel(ns);
winsize_bins = ceil(sdfac * maxjitter);
ey = compute_expectedy_lookup_table(sigma_test_list, winsize_bins);
eyssq = sum(ey .^ 2, 1);

mask = 1:numel(ns);

%identify regions in the moments where spike counts are too low to support even a single spike in the fit
pdf_sq_integral = 1 / sqrt(4 * pi * maxjitter ^ 2);
if 2 * winsize_bins + 1 <= T
    
    nssq_cs = cumsum(ns .^ 2);
    nssq_windowed = (nssq_cs(winsize_bins * 2 + 1:end) - [0; nssq_cs(1:end - winsize_bins * 2 - 1)]);
    sufficient_ns = sqrt(nssq_windowed) >= sqrt(pdf_sq_integral) / 2;
    subindex = winsize_bins + 1:T - winsize_bins;
    mask(subindex(~sufficient_ns)) = [];
    %FIXME do the edges too
    
end

%break the data into blocks to be processed seperately to reduce the dimension of optimization problems etc.
blocks = find(diff([-inf mask]) > 2 * winsize_bins);
blocke = find(diff([mask  inf]) > 2 * winsize_bins);
gfit = zeros(size(ns));
[st, stsd] = deal([]);

for blockind = 1:numel(blocks)
    
    blockmask = mask(blocks(blockind):blocke(blockind));
    block_data_index = max(1, blockmask(1) - winsize_bins):min(T, blockmask(end) + winsize_bins);
    [block_st, block_stsd, gfit(block_data_index)] = gfit_singleblock(ns(block_data_index), ey, eyssq, minjitter, maxjitter, blockmask - block_data_index(1) + 1, winsize_bins, sigma_test_list);
    st = [st; block_st + block_data_index(1) - 1];
    stsd = [stsd; block_stsd];
    
end

if ~isempty(st)
    [st, si] = sort(st);
    stsd = stsd(si);
end
warning(wstate.state, 'MATLAB:nearlySingularMatrix');


function [st, stsd, gfit] = gfit_singleblock(ns, ey, eyssq, minjitter, maxjitter, mask, winsize_bins, sigma_test_list)
T = numel(ns);
x = (0:T)' + 0.5; %bin edges between time intervals corresponding to individual counts. first count corresponds to a time of 1 bin

residual = -ns;
residual(end + 1) = 0;
e = sum(residual .^ 2);

index_mat = bsxfun(@plus, (-winsize_bins:winsize_bins)', 1:T);
index_mat(index_mat < 1 | index_mat > T) = numel(residual); %indexes the extra zero

st = []; stsd = [];

while ~isempty(mask)
    
    %try to place a new spike using a fixed list of sigma values
    [bestediff, bestsigmaind, besttind, submask] = test_sigma_vals(residual, index_mat(:, mask), ey, eyssq);
    
    if ~any(submask), break; end %couldn't reduce error function
    
    next_st = [st; mask(besttind)];
    mask = mask(submask);
    next_stsd = [stsd; sigma_test_list(bestsigmaind)];
    
    %optimize time and jitter on all spikes
    [next_st_opt, next_stsd_opt, nexte_opt, evec, iter_used] = optimize_fit(ns, next_st', next_stsd', x, minjitter, maxjitter); %FIXME mex / lbfgs!!
    
    if nexte_opt > e
        
        warning('fitting procedure for gaussian spike trains terminated earlier than expected');
        break;
        
    end
    
    residual(1:end - 1) = evec;
    st = next_st_opt;
    stsd = next_stsd_opt;
    e = nexte_opt;
    
end

gfit = residual(1:end-1) + ns;


function z = sigma2z(sigma, minjitter, maxjitter)
n = (sigma - minjitter) / (maxjitter - minjitter);
z = log(n ./ (1 - n));


function [sigma, dsigma_dz] = z2sigma(z, minjitter, maxjitter)
zexp = exp(z);
sigma = minjitter + (maxjitter - minjitter) * zexp ./ (1 + zexp);
sigma(isinf(zexp)) = maxjitter;
dsigma_dz = (sigma - minjitter) ./ (1 + zexp); %avoids NaNs


function [mu, sigma, essq, evec, iter] = optimize_fit(ns, mu0, sigma0, x, minjitter, maxjitter)

%we parameterize sigma by z, where
%sigma = minjitter + (maxjitter - minjitter) * exp(z) / (1 + exp(z))
%and
%dsigma_dz = (maxjitter - minjitter) * exp(z) / (1 + exp(z)) ^ 2

T = numel(ns);
max_iter = 20;
mintdiff = 0.01;
p = [mu0 ...
    sigma2z( min(maxjitter * 0.999, max(minjitter * 1.001, sigma0)), minjitter, maxjitter ) ...
    ]';
dp = zeros(size(p));
n = numel(mu0);
iter = 0;
elist = [];
while iter < max_iter
    
    iter = iter + 1;
    
    [evec, J] = fitpeaks_errfunc(ns, p(1:n)', p(n + 1:end)', x, minjitter, maxjitter);
    elist(1, iter) = sum(evec .^ 2);
    
    H = J' * J;
    params_used = max(abs(H)) > sqrt(eps);
    dp(:) = 0;
    dp(params_used) = H(params_used, params_used) \ (-J(:,params_used)' * evec);
    
    prev_p = p;
    prev_evec = evec;
    prev_essq = sum(evec.^2);
    
    p = prev_p + dp;
    dp = p - prev_p;
    dp(p == prev_p) = 0; %handles +/- inf
    evec = fitpeaks_errfunc(ns, p(1:n)', p(n + 1:end)', x, minjitter, maxjitter);
    essq = sum(evec .^ 2);
    pdiff = calc_pdiff(p, prev_p, n, minjitter, maxjitter);
    
    while essq > prev_essq && any(pdiff >= mintdiff)
        
        dp = dp / 2;
        p = prev_p + dp;
        evec = fitpeaks_errfunc(ns, p(1:n)', p(n + 1:end)', x, minjitter, maxjitter);
        essq = sum(evec .^ 2);
        pdiff = calc_pdiff(p, prev_p, n, minjitter, maxjitter);
        
    end
    
    if all(pdiff < mintdiff)
        
        if essq > prev_essq
            evec = prev_evec;
            essq = prev_essq;
            p = prev_p;
            %warning('failed to reduce error on final iteration'); %#ok<WNTAG>
        end
        break; %convergence succesful
        
    end
    
end

mu = p(1:n);
sigma = z2sigma(p(n + 1:end), minjitter, maxjitter);


function maxd = calc_pdiff(p, prev_p, n, minjitter, maxjitter)
maxdmu = max(abs(p(1:n) - prev_p(1:n)));
maxdsigma = max(abs( z2sigma(p(n + 1:end), minjitter, maxjitter) - z2sigma(prev_p(n + 1:end), minjitter, maxjitter) ));
maxd = max(maxdmu, maxdsigma);


function [evec, J] = fitpeaks_errfunc(ns, mu, z, x, minjitter, maxjitter) %mu, v and x are valued in time bins. mu and v should be row vectors, x a column vector
%mu and z are row vectors
[sigma, dsigma_dz] = z2sigma(z, minjitter, maxjitter);
d = bsxfun(@rdivide, bsxfun(@minus, x, mu), sigma);
P = quickncdf(d);
sP = sum(P, 2);
ens = sP(2:end) - sP(1:end-1);
evec = ens - ns;
assert(~any(isnan(evec)), 'NaN value in normal cdf difference - ns');
if nargout < 2, return; end

p = quicknpdf(d);
de_dmu = bsxfun(@rdivide, p(2:end,:) - p(1:end - 1,:), -sigma);
q = p .* bsxfun(@rdivide, d, -sigma);
de_dsigma = q(2:end, :) - q(1:end - 1, :);
de_dz = bsxfun(@times, de_dsigma, dsigma_dz); %chain rule, dsigma_dv = exp(v)
J = [de_dmu de_dz];


function [bestediff, bestsigmaind, besttind, submask] = test_sigma_vals(residual, index_mat, ey, eyssq)
resvalmat = residual(index_mat);
ediff = -bsxfun(@plus, 2 * ey' * resvalmat, eyssq'); %difference in e^2 after - before adding new component
submask = any(ediff > 0, 1);
[bestediff_allsigma, besttind_allsigma] = max(ediff, [], 2);
[bestediff, bestsigmaind] = max(bestediff_allsigma);
besttind = besttind_allsigma(bestsigmaind);


function ey = compute_expectedy_lookup_table(sigma_test_list, winsize_bins)
ey = nan(2 * winsize_bins + 1, numel(sigma_test_list));
for u = 1:numel(sigma_test_list)
    sigma = sigma_test_list(u);
    xvals = ((-(winsize_bins + 1):winsize_bins) + 0.5) / sigma;
    ncdf = quickncdf(xvals);
    ey(:, u) = ncdf(2:end) - ncdf(1:end - 1);
end