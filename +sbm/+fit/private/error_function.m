function [e, grad, fhat, states, eseg, R, states_eq, eneuron] = error_function( ...
    Pfunc, simulation_function, ...
    f, fBL, spikecounts, stepsize_eachsubseg, nsteps_eachsubseg, ...
    n_extrabins_eachsubseg, p, pmin, isobs, symgrad, ...
    segmentind, p_seglist, gradreldiff, gradmindiff, neuronindex_eachsubseg, focalplaneindex_eachsubseg, paramopts, states, states_eq, redostates, subseglist)
%this is the error function. it also calculates its own gradients numerically

P_eachseg = Pfunc(p);
P_eachsubseg = P_eachseg(segmentind);

nsubseg = numel(f);

if redostates

    %simulate binding/unbinding based on spike counts
    [states, states_eq] = simulation_function( ...
        stepsize_eachsubseg, spikecounts, P_eachsubseg, states, subseglist);

    assert(~any(cellfun(@(v) any(isnan(v(:))), states)), 'simulation failed');

end

[e, eseg, R, fhat, eneuron] = fitR(numel(subseglist), states(subseglist), f(subseglist), fBL(subseglist), states_eq(:, subseglist), P_eachsubseg(subseglist), ...
    n_extrabins_eachsubseg(subseglist), nsteps_eachsubseg(subseglist), paramopts, neuronindex_eachsubseg(subseglist), focalplaneindex_eachsubseg(subseglist));

if any(isnan(eseg))

    R = []; e = inf; grad = nan(numel(p), 1); fhat = cell(1, nsubseg);
    return;

end

if nargout < 2

    return;

end

%calculate gradients
[e_plus, e_plus_plus, e_minus, e0, grad] = deal(nan(numel(p), 1));
for j = 1:numel(p)

    %numeric gradient. note fdc changes automatically as we perturb other parameters
    subseglist_forgrad = find(ismember(segmentind, p_seglist{j}));
    if any(focalplaneindex_eachsubseg(subseglist_forgrad)) %need to optimize over all segements that share an fpi

        fpilist = unique(focalplaneindex_eachsubseg(subseglist_forgrad));
        fpilist(fpilist == 0) = [];
        subseglist_forgrad = unique([subseglist_forgrad find(ismember(focalplaneindex_eachsubseg(:)', fpilist))]);

    end

    %recalculate states only if necessary
    redostates_forgrad = ~isobs(j);

    %fixme check pmax?
    
    %perturb element j of the parameter vector
    delta = max(gradmindiff, abs(p(j)) * gradreldiff);
    p_plus = p;
    p_plus(j) = p_plus(j) + delta;
    
    %copy the variables so that the originals don't get modified by mex files
    states_copy = copy_matlab_var(states);
    states_eq_copy = copy_matlab_var(states_eq);    
    
    e_plus(j) = error_function(...
        Pfunc, simulation_function, ...
        f, fBL, spikecounts, stepsize_eachsubseg, nsteps_eachsubseg, ...
        n_extrabins_eachsubseg, p_plus, pmin, isobs, symgrad, ...
        segmentind,  ...
        [], gradreldiff, gradmindiff, neuronindex_eachsubseg, focalplaneindex_eachsubseg, paramopts, states_copy, states_eq_copy, redostates_forgrad, subseglist_forgrad);
    
    if symgrad
        
        % we don't need a 2nd copy of states/states_eq since, everything
        % that was modified in the positive delta direction should be
        % modified again. when the gradient is later computed w.r.t. other
        % parameters, the necessary segments will also be re-run
        
        if p(j) - delta > pmin(j)
            
            p_minus = p;
            p_minus(j) = p(j) - delta;
            e_minus(j) = error_function(...
                Pfunc, simulation_function, ...
                f, fBL, spikecounts, stepsize_eachsubseg, nsteps_eachsubseg, ...
                n_extrabins_eachsubseg, p_minus, pmin, isobs, symgrad, ...
                segmentind,  ...
                [], gradreldiff, gradmindiff, neuronindex_eachsubseg, focalplaneindex_eachsubseg, paramopts, states_copy, states_eq_copy, redostates_forgrad, subseglist_forgrad);
            
            grad(j) = (e_plus(j) - e_minus(j)) / (p_plus(j) - p_minus(j));
            
        else
            
            p_plus_plus = p;
            p_plus_plus(j) = p(j) + 2 * delta;
            e_plus_plus(j) = error_function(...
                Pfunc, simulation_function, ...
                f, fBL, spikecounts, stepsize_eachsubseg, nsteps_eachsubseg, ...
                n_extrabins_eachsubseg, p_plus_plus, pmin, isobs, symgrad, ...
                segmentind,  ...
                [], gradreldiff, gradmindiff, neuronindex_eachsubseg, focalplaneindex_eachsubseg, paramopts, states_copy, states_eq_copy, redostates_forgrad, subseglist_forgrad);
            
            % calculate previous error for relevant neurons only
            ntimepoints_eachsubseg = cellfun(@numel, f(subseglist_forgrad));
            e0(j) = calc_error(eseg(subseglist_forgrad), ntimepoints_eachsubseg, neuronindex_eachsubseg(subseglist_forgrad));
            
            % get slope from tangent to quadratic
            grad(j) = (-3 * e0(j) + 4 * e_plus(j) - e_plus_plus(j)) / (2 * delta);  % slope at x0 from 3-point quadratic fit
            
        end
        
    else

        % calculate previous error for relevant neurons only
        ntimepoints_eachsubseg = cellfun(@numel, f(subseglist_forgrad));
        e0(j) = calc_error(eseg(subseglist_forgrad), ntimepoints_eachsubseg, neuronindex_eachsubseg(subseglist_forgrad));
        
        grad(j) = (e_plus(j) - e0(j)) / (p_plus(j) - p(j));

    end
    
end


function [e, eseg, R, fhat, eneuron] = fitR(nsubseg, states, f, fBL, states_eq, P_eachsubseg, ...
    n_extrabins_eachsubseg, nsteps_eachsubseg, paramopts, neuronindex_eachsubseg, focalplaneindex_eachsubseg)
%calculate relative fluorescence changes over time and at equilibrium. discard states one full imaging time step or more before the first fluorescence measurement
%this function fitsR or uses the current param values as appropriate based on paramopts.autoR
if ismember('R', paramopts.focalplanespecific)

    Rind = focalplaneindex_eachsubseg;

elseif ismember('R', paramopts.neuronspecific)

    Rind = neuronindex_eachsubseg;
    
else
    
    Rind = ones(nsubseg, 1);

end
nfp = max(Rind(:));
ntimepoints_eachsubseg = cellfun(@numel, f);
maxdff = cell(1, nsubseg); %value of (F / F_baseline) when R = 0
brightnessincrease_eq = nan(1, nsubseg);
for u = 1:nsubseg

    S = P_eachsubseg(u).S;
    states_at_fobs = states{u}(:, n_extrabins_eachsubseg(u) + nsteps_eachsubseg(u):nsteps_eachsubseg(u):end);
    
    db = P_eachsubseg(u).dbrightness;
    ii = 3:3 + numel(db) - 1;
    brightnessincrease = db * states_at_fobs(ii, :) / S;  % first state is calcium, second is ground state
    brightnessincrease_eq(u) = db * states_eq(ii, u) / S;  % first state is calcium, second is ground state
    maxdff{u} = (brightnessincrease - brightnessincrease_eq(u))' / (brightnessincrease_eq(u) + 1); %fractional brightness increase over baseline

end

brightnessincrease_eq(isnan(brightnessincrease_eq)) = 0;  %e.g. S = 0

[R, v] = deal(nan(nfp, 1));
for jj = 1:nfp

    ii = Rind == jj;
    if ~any(ii), continue; end

    if ~isfield(P_eachsubseg, 'dbrightness')
        
        beq = brightnessincrease_eq(find(ii, 1));
        assert(all(brightnessincrease_eq(ii) == beq), 'equilibrium brightness must be the same for all segments with the same focal plane');
        
    end

    next_f   = cat(1,      f{ii});
    next_fBL = cat(1,    fBL{ii});
    next_dff = cat(1, maxdff{ii});

    if ~paramopts.autoR

        R(jj) = P_eachsubseg(find(ii, 1)).R;        

    elseif any(isnan(next_dff)) || all(next_dff == 0) %this can happen if e.g. the indicator never gets bright for these model parameters since some on rate is zero

        [v(jj), R(jj)] = deal(nan);
        continue;

    elseif paramopts.maxRval == paramopts.minRval

        R(jj) = paramopts.minRval;        

    else

        v(jj) = (next_fBL .* next_dff) \ (next_f - next_fBL);
        
        if isfield(P_eachsubseg, 'dbrightness')
            
            R(jj) = min(paramopts.maxRval, max(paramopts.minRval, (1 - v(jj)) / v(jj)));
            
        else
            
            R(jj) = min(paramopts.maxRval, max(paramopts.minRval, beq * (1 - v(jj)) / v(jj)));
            
        end
        

    end
    
    if isfield(P_eachsubseg, 'dbrightness')
        
        v(jj) = 1 / (1 + R(jj));
        
    else
        
        v(jj) = beq / (beq + R(jj));
        
    end

end

%calculated model-predicted fluorescence and error
fhat = cell(1, nsubseg);
eseg =  nan(1, nsubseg);
for u = 1:nsubseg

    if isnan(v(Rind(u))) %this can happen if e.g. the indicator never gets bright for these model parameters since some on rate is zero

        fhat{u} = fBL{u};

    else

        fhat{u} = fBL{u} .* (maxdff{u} * v(Rind(u)) + 1);

    end
    eseg(u) = sum((f{u} - fhat{u}) .^2);

end
[e, eneuron] = calc_error(eseg, ntimepoints_eachsubseg, neuronindex_eachsubseg);


function [e, eneuron] = calc_error(eseg, ntimepoints, neuronindex)
nlist = unique(reshape(neuronindex, 1, []));
eneuron = zeros(1, numel(nlist));
for k = 1:numel(nlist)
    
    n = nlist(k);
    ii = neuronindex == n;
    eneuron(k) = sum(eseg(ii)) / sum(ntimepoints(ii));    

end
e = sum(eneuron);