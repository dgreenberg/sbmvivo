function [M, P, V, lik, prevP, x, moments, prevlik] = run_alg(varargin)
% [M, P, V, lik, prevP, x, moments, prevlik] = run_alg(F,dt,indicator,V,P,x,n,nA2D,st,settingsindex,sessionindex,focalplaneindex)
% this function runs sequential monte carlo spike inference with parameter detection on a fluorescence time-series
% and outputs the inferred distributions and parameter estimates
% Optional inputs can be either omitted or left empty.
%
% Inputs
% F: fluorescence time series. row vector or cell array of row vectors
% dt: time step for F (sec)
% indicator: a string indicating the calcium indicator used to image this neuron
% V: structure of options
% P: structure of initial parameter estimates
% x: stimulus for x (optional), with numel(F) colums. one row should consist of ones to allow stimulus independent firing
% n: spikes per image frame (optional), will be used for initial param est. when n is nonempty and V.use_true_n == 1, posterior on hidden states is always calculated given this spike train
% nA2d: number of analog to digital measurements per fluorescence value
%
% Outputs
% M: structure containing mean, variance, and percentiles of inferred distributions
% P: structure array containing parameter estimates
% V: structure Variables for algorithm to run

%FIMXE -- document all inputs and outputs
[F,V0,P0,dt,indicatorstring,x,true_n,nA2D,st] = parse_alg_inputs(varargin);
[P,V] = sbm.init.initialize_params_and_settings_from_data(F, dt, P0, V0, indicatorstring, nA2D); %assigns any parameters and settings not supplied by the user using default values, some of which may be data-dependent
prevP = P([],:); prevlik = repmat(empty_likstruct, 0, numel(F));

do_paramest = V.min_iter_paramest > 0 && ~strcmpi(V.parameterestimation, 'none');
if do_paramest
    
    [P, prevP, prevlik] = estimate_parameters(V, P, dt, F, nA2D);
    [moments, lik, ~] = filter_and_smooth(F, dt, V, P, nA2D, st);
    
elseif V.profile_firing_rate || V.profile_zeta

    prevP = [prevP; P];
    %if we're profiling but not estimating params, profile over FR now.
    [~, lik, ~, P, moments] = profile_likelihood(F, dt, V, P, nA2D, st);
    
else
    
    [moments, lik, ~, ~] = filter_and_smooth(F, dt, V, P, nA2D, st);
    
end
if V.gfitst
    
    M = gfitspiketrain(moments, V, dt);
    
else
    
    M = struct();
    
end
    
if (do_paramest || V.profile_firing_rate || V.profile_zeta) && V.verbose > 2
    
    fprintf('\nFinal parameter changes:\n');
    print_param_changes([prevP(1,:); P], V);
    
    if ~isempty(prevlik)
        
        print_likelihood(prevlik(1,:), 'LL at start  ');
        
    end
    
    print_likelihood(lik,          'LL at end    ');
    
end


function [F,V,P,dt,indicatorstring,x,true_n,nA2D,st] = parse_alg_inputs(invars)
nvars = numel(invars);
assert(nvars > 0,'fluorescence must be provided as the first input');
F = invars{1};
if isnumeric(F)
    
    F = {F};
    
else
    
    assert(iscell(F),'F must be a numeric row vector or cell array of numeric row vectors');
    F = reshape(F, 1, []);
    
end
assert(nvars > 1, 'time step must be provided');
dt = invars{2};
assert(isnumeric(dt) && size(dt,1) == 1 && isvector(dt) && all(dt > 0),'dt must be a scalar or row vector and greater than 0');
dt = double(dt);
if numel(dt) == 1
    dt = repmat(dt, 1, numel(F));
end
if nvars < 3 || isempty(invars{3})
    indicatorstring = '';
else
    indicatorstring = invars{3};
end
assert(isa(indicatorstring,'char'),'indicator input must be a string');
if nvars < 4 || isempty(invars{4})
    V = struct;
else
    V = invars{4};
end
assert(numel(V) == 1, 'V is incorrectly sized');
if nvars < 5 || isempty(invars{5})
    P = struct;
else
    P = invars{5};
end
if numel(P) == 1
    P = repmat(P, 1, numel(F));
else
    assert(numel(P) == numel(F),'P is incorrectly sized');
    P = reshape(P, 1, numel(F));
end
if nvars < 6 || isempty(invars{6})
    x = cell(1, numel(F));
else
    x = invars{6};
    if isnumeric(x)
        x = {x};
    else
        assert(iscell(x) && numel(x) == numel(F));
        x = reshape(x, 1, []);
    end
end
if nvars < 7 || isempty(invars{7})
    true_n = cell(1, numel(F));
else
    true_n = invars{7};
    if isnumeric(true_n)
        true_n = {true_n};
    else
        true_n = reshape(true_n, 1, []);
    end
    assert(iscell(true_n) && numel(true_n) == numel(F) && size(true_n,1) == 1 && size(true_n,2) == numel(F), 'invalid input for n');
end
if nvars < 8 || isempty(invars{8})
    nA2D = ones(1, numel(F));
else
    nA2D = checkvectorinput(invars{8}, numel(F), 'nA2D');
    if all(isnan(nA2D))
        nA2D(:) = 1;
    elseif any(isnan(nA2D))
        warning('sbm:nA2DisNaN','NaN values for nA2D have been replaced with the mean of the non-nan values');
        nA2D(isnan(nA2D)) = mean(nA2D(~isnan(nA2D)));
    end
end
if nvars < 9 || isempty(invars{9})
    st = cell(1, numel(F));
else
    st = invars{9};
    for u = 1:numel(F)
        if ~isempty(st{u})
            assert(all(isnumeric(st{u}) & isfinite(st{u}) & ~any(isnan(st{u})) & ~any(imag(st{u}))),'invalid spike time input');
        end
    end
end
assert(nvars < 10,'maximium of 9 inputs');

for u = 1:numel(F) %check that each element of each cell array is correct
    
    assert(isnumeric(F{u}) && isvector(F{u}) && ~isempty(F{u}), 'F must be a numeric vector or cell array of numeric vectors');
    F{u} = double(reshape(F{u}, 1, []));        
    T = numel(F{u}); %number of observations
    assert(~any(isnan(F{u})) && isreal(F{u}),'nan and complex values for F are not supported')
    
    if isempty(x{u})
        x{u} = ones(1, T);
    end
    assert(numel(size(x{u})) == 2 && size(x{u},2) == T,'x and F do not match');
    assert(~any(isnan(x{u})) && isreal(x{u}),'nan and complex values for x are not supported');
    
    if ~isempty(true_n{u})
        assert(numel(size(true_n{u})) == 2 && size(true_n{u},2) == T && size(true_n{u},1) == 1,'true_n and F do not match');
        assert(~any(isnan(true_n{u})) && isreal(true_n{u}),'nan and complex values for true_n are not supported');
    end
    
end
xdim = cellfun(@(vv) size(vv,1), x);
assert(all(xdim == xdim(1)), 'the dimension of x must be the same for all segments so that a single filter can be learned');