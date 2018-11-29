function [b, tau_s, s, blmask, beta] = fit_unbinding_dynamics_withoffset(f, it, st, tau_s0, lambda, tafterc0, tafters0, tbefore, minbeta, minblmaskpoints)
%[b, tau_s, s, blmask] = fit_unbinding_dynamics(f, it, st, tau_s0, lambda, tafterc0, tafters0, tbefore)
%this function fits the dye time constant by examing periods where free calcium has fallen to near zero but unbinding has not yet been completed
%it can also incorporate non-baseline fluorescence measurements into accurate estimates of the baseline value toward which decay is occurring
%an offset is inferred. b is never negative

%the fact we put a maximum time constant on the decay is critical, as otherwise we would be able to get near-zero error by using a huge s and a small, arbitrary b

%parse inputs
f = reshape(f, [], 1);
st = reshape(st, 1, []);
dt = mean(diff(it));
if ~exist('lambda','var') || isempty(lambda)
    
    t_drift_1sd = 240; %time in seconds for baseline drift to cause a fluorsecence change equal to the data median
    lambda = t_drift_1sd / dt;
    
end
if ~exist('tbefore','var')  || isempty(tbefore),  tbefore = 1.5 * dt; end
if ~exist('tafters0','var') || isempty(tafters0), tafters0 = 6;       end
if ~exist('tafterc0','var') || isempty(tafterc0), tafterc0 = 5.5;     end
if ~exist('tau_s0','var')   || isempty(tau_s0),   tau_s0 = 0.9;       end
if ~exist('minbeta','var')  || isempty(minbeta),  minbeta = -inf;     end
if ~exist('minblmaskpoints', 'var') || isempty(minblmaskpoints), minblmaskpoints = 1; end

%get an initial estimate of baseline
[blmask_init, blmask] = deal(true(size(f))); %true vector in the size of f in both variables
for t = st
    
    blmask_init(it >= t - tbefore & it <= t + tafters0) = false; %set false in tbefore-tafter0 environment of spike times
    blmask     (it >= t - tbefore & it <= t + tafterc0) = false;
    
end
if mean(blmask_init) < 0.05 %a lot of spikes?
    
    b = medfilt_windowed(f, round(60 / dt)); %smooth f over a 30s region, median helps to preserve edges?
    
else
    
    gwslow = 5 / dt;
    gwfast = 0.5 / dt;
    b = (gauss_filter(f .* double(blmask_init), gwslow) + gauss_filter(f .* double(blmask_init), gwfast)) ...
        ./ (gauss_filter(double(blmask_init), gwslow) + gauss_filter(double(blmask_init), gwfast)); %Gauss filter
    
end
borig = b;
blmask = blmask & ~isnan(b);
assert(sum(blmask) >= minblmaskpoints, 'can''t find baseline segment'); %error when everything is signal range
if ~any(blmask)
    
    tau_s = nan;
    s = nan(size(b));
    beta = 0;
    return;
    
end

%divide data into segments
segments = [find([true; ~blmask(1:end-1)] & blmask) ...
    find(blmask & [~blmask(2:end); true])]; %find indexes of edges of tbefore-tafter0 environment of spike times (in pairs start/begining of zones between environments)

a = exp(-dt / tau_s0);          %convert initial estimate of decay constant to multiplicative factor per time step
amin = min(a, exp(-dt / 0.3));  %minimal and maximal decay rate for optimisation
amax = max(a, exp(-dt / 1.5));
beta = min(0, min(f)); %initialize offset
b = b - beta; %should make b nonnegative

[Q, D] = difssq_penaltymat(segments);

fmasked = f(blmask); %part of flourescence timeseries outside environments (where c is supposed to be decayed)
b(~blmask) = nan;

if sum(double(~isnan(b))) < 3
    
    s = b * 0; tau_s = 1;
    warning('non enough data to accurately estimate decay constant etc.');
    return;
    
end

ws = warning('query', 'optim:quadprog:WillRunDiffAlg');
warning('off', 'optim:quadprog:WillRunDiffAlg');
e = inf; efac = 1e-5; iter = 0; elist = [];
while true %we will do 3-way coordinate descent on b, beta and the decay factor a
    
    eprev = e;
    [~, s, a] = fitdecayrate(f - beta, b, segments, amin, amax, a); %don't return error here since it's only for the data mismatch term, not the smoothness term
    beta = min(min(f), max(minbeta, mean_nonnan(f - b .* (1 + s))));
    if isempty(elist), sorig = s; end
    [e, b(blmask)] = fitbl(fmasked - beta, s(blmask), lambda, Q, D);
    beta = min(min(f), max(minbeta, mean_nonnan(fmasked - b(blmask) .* (1 + s(blmask)))));
    e = sum((fmasked - b(blmask) .* (1 + s(blmask)) - beta).^2) + lambda * b(blmask)' * Q * b(blmask);
    elist(end+1) = e; %#ok<AGROW>
    iter = iter + 1;
    if e < eps || (eprev - e) / eprev < efac
        break;
    end
    
end
warning(ws.state, 'optim:quadprog:WillRunDiffAlg');
tau_s = -dt / log(a);
b = blgapfill(b, segments);


function [Q, D] = difssq_penaltymat(segments)
segpts = diff(segments,1,2)+1; %elements/points of zones between environments
Q = sparse([], [], [], sum(segpts), sum(segpts), sum(1 + 3 * (segpts - 1))); %sum(segpts)xsum(segpts) zero sparse matrix with space allocated for [last parameter]
D = sparse([], [], [], sum(segpts), sum(segpts), 2 * (sum(segpts) - 1));
nprev = 0;
for k = 1:size(segments,1) %for each segment
    nextn = segments(k,2) - segments(k,1) + 1; %==segpts(k)?
    for j = 1:nextn - 1
        
        Q(nprev + j, nprev + j) = Q(nprev + j, nprev + j) + 1; %#ok<*SPRIX>
        Q(nprev + j + 1, nprev + j + 1) = Q(nprev + j + 1, nprev + j + 1) + 1;
        Q(nprev + j, nprev + j + 1) = -1;
        Q(nprev + j + 1, nprev + j) = -1;
        
        D(nprev + j, nprev + j:nprev + j + 1) = [-1 1];
        
    end
    if k < size(segments, 1) %every but the last segment
        
        nbinsdiff = segments(k+1,1) - segments(k,2);
        Q(nprev + nextn, nprev + nextn) = Q(nprev + nextn, nprev + nextn) + 1 / nbinsdiff; %#ok<*SPRIX>
        Q(nprev + nextn + 1, nprev + nextn + 1) = Q(nprev + nextn + 1, nprev + nextn + 1) + 1 / nbinsdiff;
        Q(nprev + nextn, nprev + nextn + 1) = -1 / nbinsdiff;
        Q(nprev + nextn + 1, nprev + nextn) = -1  / nbinsdiff;
        
        D(nprev + nextn, nprev + nextn:nprev + nextn + 1) = [-1 1] / sqrt(nbinsdiff);
        
    end
    nprev = nprev + nextn;
end


function b = blgapfill(b, segments)
b(1:segments(1,1)-1) = b(segments(1,1));
for k = 1:size(segments,1) - 1
    ii = segments(k,2) + 1:segments(k+1,1) - 1;
    b(ii) = b(segments(k,2)) + (b(segments(k+1,1)) - b(segments(k,2))) * (1:numel(ii)) / (numel(ii) + 1);
end
b(segments(end,2)+1:end) = b(segments(end,2));


function [e, bmasked] = fitbl(fmasked, smasked, lambda, Q, D)
L = speye(numel(smasked));
L(1:size(L,1) + 1:end) = (1 + smasked) .^ 2;
A = lambda * Q + L;
d = double(fmasked .* (1 + smasked));
bmasked = A \ d; %solve by setting derivative of error function to zero
if any(bmasked < -eps) %need to solve a constrained quadratic program. this should happen rarely, if ever.
    
    bmasked = quadprog(2 * A, -2 * d, [], [], [], [], zeros(numel(fmasked), 1), [], [], optimset('display','off'));
    
else
    
    bmasked = max(bmasked, 0);
    
end
e = sum((fmasked - bmasked .* (1 + smasked)).^2) + lambda * bmasked' * Q * bmasked;


function [e, s, a] = fitdecayrate(f, b, segments, amin, amax, aprev)
ofunc = @(a) fittransients(f, b, segments, a);
%a = fminbnd(ofunc,amin,amax,optimset('display','off','tolx',1e-5));
[a,fv,ef,op] = fminbnd(ofunc,amin,amax);
[e, s] = fittransients(f, b, segments, a);
eprev = ofunc(aprev);
if e > eprev
    %warning('failed to decrease error, using previous decay rate'); %#ok<WNTAG>
    a = aprev;
    [e, s] = fittransients(f, b, segments, a);
end


function [e, s] = fittransients(f, b, segments, a)
e = 0;
s = nan(size(f));
for k = 1:size(segments, 1) % for every segment
    ii = (segments(k,1):segments(k,2))';
    apowers = a .^ (ii - ii(1));
    decaytrans = b(ii) .* apowers;
    smax = max(0, decaytrans \ (f(ii) - b(ii))); %least squares solutions for decaytrans*x=f(ii) - b(ii)
    s(ii) = apowers * smax;
    e = e + sum((f(ii) - b(ii) .* (1 + s(ii))).^2);
end