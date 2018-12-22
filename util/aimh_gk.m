function [samples, targetlogpvals, acceptance_ratio, accepted, gstar, extradata, first_update, preliminary_phase_end] = aimh_gk( ...
    targetlogpfunc, proposal, p0, temperature, mindraws, maxdraws, halt_on_n_accepted, displayfunc, box_constraints, w1, w2, L, alpha_L, M, alpha_M, w1_init)
%proposal must be a gmix struct
%
%http://amstat.tandfonline.com/doi/abs/10.1198/jcgs.2009.07174
update_every = 1;
k1 = 9; %25
k2 = 9; %16

if ~exist('p0', 'var') || isempty(p0)
    samples = gmix_sample(proposal);
else
    samples = p0(:)';
end
datadim = size(samples, 2);
[targetlogpvals, extradata{1}] = targetlogpfunc(samples(1, :));

if ~exist('temperature', 'var') || isempty(temperature)
    temperature = 1;
end
if ~exist('w1', 'var') || isempty(w1)
    w1 = 0.05;
end
if ~exist('w2', 'var') || isempty(w2)
    w2 = 0.15;
end
if ~exist('L', 'var') || isempty(L)
    L = 10;
end
if ~exist('alpha_L', 'var') || isempty(alpha_L)
    alpha_L = 0.1;
end
if ~exist('M', 'var') || isempty(M)
    M = 20;
end
if ~exist('alpha_M', 'var') || isempty(alpha_M)
    alpha_M = 0.02;
end
if ~exist('w1_init', 'var') || isempty(w1_init)
    w1_init = 0.75; %0.6
end
if ~exist('mindraws', 'var') || isempty(mindraws)
    mindraws = 400;
end
if ~exist('maxdraws', 'var') || isempty(maxdraws)
    maxdraws = mindraws * 3;
end
if ~exist('halt_on_n_accepted', 'var') || isempty(halt_on_n_accepted)
    halt_on_n_accepted = max(10, ceil(mindraws / 4));
end
if ~exist('displayfunc', 'var') || isempty(displayfunc)
    displayfunc = [];
end
if ~exist('box_constraints', 'var') || isempty(box_constraints)
    box_constraints = repmat([-inf; inf], 1, numel(p0));
end

[first_update, preliminary_phase_end] = deal([]);

min_accepted_draws_forupdate = 5 * datadim;
max_accepted_draws_preliminary_phase = 20 * datadim;

g0 = mix_mixtures(proposal, widen_tails(proposal, k1), w1_init);
q = g0;
gstar = [];

[acceptance_ratio, accepted] = deal([]);
[ndraws, naccepted] = deal(0);
preliminary_phase = true;
t0 = now;

preliminary_phase_completed = false;
while ndraws < maxdraws && (ndraws < mindraws || sum(accepted) < halt_on_n_accepted)
    
    next_draw = gmix_sample(q);
    while any(next_draw < box_constraints(1, :)) || any(next_draw > box_constraints(2, :))
        next_draw = gmix_sample(q);
    end    
    
    [nexttargetlogpval, next_extradata] = targetlogpfunc(next_draw);
    
    logdiff = nexttargetlogpval - targetlogpvals(end) + gmix_logpdf(q, samples(end, :)) - gmix_logpdf(q, next_draw);
    acceptance_ratio(1, end + 1) = min(1, exp(logdiff / temperature));
    accepted(1, end + 1) = rand < acceptance_ratio(1, end);
    
    if accepted(1, end)
        
        naccepted = naccepted + 1;
        samples(end + 1, :) = next_draw;
        targetlogpvals(1, end + 1) = nexttargetlogpval;
        extradata{1, end + 1} = next_extradata;
        
    else
        
        samples(end + 1, :) = samples(end, :);
        targetlogpvals(1, end + 1) = targetlogpvals(1, end);
        extradata{1, end + 1} = extradata{1, end};
        
    end
    ndraws = ndraws + 1;
    
    update_mixture = (~mod(ndraws, update_every) || ndraws == maxdraws || ndraws == mindraws) && naccepted >= min_accepted_draws_forupdate;
    
    if preliminary_phase && naccepted >= min_accepted_draws_forupdate
        
        preliminary_phase_completed = ...
            (ndraws >= M && all(acceptance_ratio(end - M + 1:end) >= alpha_M)) || ...
            naccepted >= max_accepted_draws_preliminary_phase;
        
        if naccepted == min_accepted_draws_forupdate || ...
                preliminary_phase_completed || ...
                (ndraws >= L && mean(accepted(end - L + 1:end)) < alpha_L) %too many rejections
            
            update_mixture = true;
            
        end
        
    end
    
    
    if update_mixture %update the mixture approximation to the accepted draws we have so far
        
        if isempty(gstar)
            first_update = numel(accepted);
        end
        
        if ~isempty(first_update) && size(samples, 1) - first_update >= 100 && sum(accepted(first_update + 1:end)) > 20
        
            gstar = gmfit_gk(samples(first_update + 1:end, :));
            
        else
            
            gstar = gmfit_gk(samples);
            
        end
        
        gstar = gmfit_gk(samples);
        
        if preliminary_phase && preliminary_phase_completed
            
            preliminary_phase = false;
            g0 = mix_mixtures(gstar, widen_tails(gstar, k1), w1_init);
            preliminary_phase_end = numel(accepted);
            
        end
        
        %calculate a mixture of the fit and a widened version
        g = mix_mixtures(widen_tails(gstar, k2), gstar, w2 / (1 - w1)); %based on eq. 2.7. this seems to be incompatible with the definition of gn in the next section however
        %calculate our sampling distribution:
        q = mix_mixtures(g0, g, w1);
        
    end
    
    tvec = datevec(now - t0);
    totalt_sec = tvec(end-2:end) * [3600 60 1]';
    
    if ~isempty(displayfunc)
        
        displayfunc(samples, targetlogpvals, acceptance_ratio, accepted, gstar, next_draw, extradata, preliminary_phase, totalt_sec, maxdraws);
        
    end
    
end
if ~preliminary_phase_completed
   
    warning('preliminary sampling phase did not complete');
    
end
if isempty(gstar)
    
    gstar = gmfit_gk(samples);
    warning('TM may have produced an insufficent number of samples');
    
end


function LL = gmix_logpdf(gmix, y)
%calculate log-likelihood of a single data point given this Gaussian mixture model
nclusters = numel(gmix.p);
datadim = size(y, 2);
LL_eachcomponent = nan(nclusters, 1);

for k = 1:nclusters
    
    d = y - gmix.m(k, :);
    dtCinvd = d * (gmix.V(:, :, k) \ (d'));
    LL_eachcomponent(k) = -0.5 * (dtCinvd + log(det(gmix.V(:, :, k)) * (2 * pi) ^ datadim)) + log(gmix.p(k));
    
end

LL_max = max(LL_eachcomponent); %LL from best cluster center for each data point

LL = log(sum(exp(LL_eachcomponent - LL_max))) + LL_max;


function y = gmix_sample(gmix)
i = find(rand <= cumsum(gmix.p), 1);
if isempty(i), i = numel(gmix.p); end %numerical error etc.
x = randn(size(gmix.m, 2), 1);
y = (chol(gmix.V(:,:, i)) * x)' + gmix.m(i, :);


function gmix = widen_tails(gmix, k)
gmix.V = gmix.V * k;


function gmix = mix_mixtures(gmix, gmix2, w1)
gmix.m = [gmix.m; gmix2.m];
gmix.V = cat(3, gmix.V, gmix2.V);
gmix.p = [gmix.p * w1, gmix2.p * (1 - w1)];