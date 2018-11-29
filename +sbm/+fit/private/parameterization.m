function [p0, pi, p_seglist_states, pmin, pmax, isobs] = parameterization( ...
    paramopts, P, parameter_names, nseg, neuronindex, sessionindex, focalplaneindex)

neuronindex = reshape(neuronindex, [], 1);
[p0, pmin, pmax] = deal([]);
isobs = false(0,1);

pi = struct();  %structure of parameter index vectors into p/L for each parameter name. each vector has nseg elements
p_seglist_states = {}; %list of segments which for which binding states depend on a given element of p
ii = 0;
for j = 1:numel(parameter_names)
    if ~isfield(P(1), parameter_names{j})
        fprintf('field not found, skipping: %s\n', parameter_names{j});
        continue;
    end
    
    nelements = numel(P(1).(parameter_names{j})); %how many values per instance of this parameter
    
    if nelements == 0 %parameter not used
        
        groupindex = zeros(size(neuronindex));
    
    elseif ismember(parameter_names{j}, paramopts.sessionspecific)
        
        groupindex = sessionindex;
        
    elseif ismember(parameter_names{j}, paramopts.focalplanespecific)
        
        groupindex = focalplaneindex;
        
    elseif ismember(parameter_names{j}, paramopts.neuronspecific)
        
        groupindex = neuronindex;
        
    else %global parameter
        
        groupindex = ones(nseg, 1);
        
    end
    
    ngroups = max(groupindex);
    
    %extract parameter values to go in parameter vector
    pv = nan(ngroups, nelements);
    for k = 1:ngroups
        
        pv(k, :) = P(find(groupindex == k, 1)).(parameter_names{j});
        
    end
    
    %note that when using logarithms, a positive parameter initialized at zero will be forever fixed
    mindiff_init = paramopts.mindiff_init; %default value
    [nextpmin, nextpmax] = deal(-inf, inf);
    if ismember(parameter_names{j}, {'Btot'})
        
        nextpmin = 0;
        nextpmax = inf;
        
    elseif ismember(parameter_names{j}, {'tau_B'})
        
        nextpmin = paramopts.min_tau_B;
        nextpmax = paramopts.max_tau_B;
        
    elseif ismember(parameter_names{j}, {'tau_ex'})
        
        nextpmin = paramopts.min_tau_ex;
        nextpmax = paramopts.max_tau_ex;
        
    elseif ismember(parameter_names{j}, {'cafree_normalized'})
        
        nextpmin = 0.0;
        nextpmax = 1.0;
        
    elseif ismember(parameter_names{j}, {'fratio'})
        
        nextpmin = paramopts.minfratio;
        nextpmax = paramopts.maxfratio;
        
    elseif ismember(parameter_names{j}, {'brightness_ratio'})
        
        nextpmin = 0.0;
        nextpmax = 1.0;
        
    elseif ismember(parameter_names{j}, {'k50_first'})
        
        nextpmin = paramopts.mink50_first;
        nextpmax = paramopts.maxk50_first;
        
    elseif ismember(parameter_names{j}, {'dk50'})
        
        nextpmin = paramopts.min_dk50;
        nextpmax = paramopts.max_dk50;
        
    elseif ismember(parameter_names{j}, {'tau_mm'})
        
        nextpmin = paramopts.min_tau_mm;
        nextpmax = paramopts.max_tau_mm;
    
    elseif ismember(parameter_names{j}, {'koff'})
        
        nextpmin = 1.0 / paramopts.maxtauval;
        nextpmax = 1.0 / paramopts.mintauval;
        if ~isinf(paramopts.max_rate_constant_adjustment)
            
            nextpmin = max(nextpmin, pv * (1 - paramopts.max_rate_constant_adjustment));
            nextpmax = min(nextpmax, pv * (1 + paramopts.max_rate_constant_adjustment));
            
        end
        
    elseif ismember(parameter_names{j}, {'kon'})
        
        nextpmin = paramopts.minkon;
        nextpmax = paramopts.maxkon;        
        
        if ~isinf(paramopts.max_rate_constant_adjustment)
            
            nextpmin = max(nextpmin, pv * (1 - paramopts.max_rate_constant_adjustment));
            nextpmax = min(nextpmax, pv * (1 + paramopts.max_rate_constant_adjustment));
            
        end
        
    elseif ismember(parameter_names{j}, paramopts.KDvals)
        
        nextpmin = paramopts.minKDval;
        nextpmax = paramopts.maxKDval;        
        
    elseif ismember(parameter_names{j}, {'A'})
        
        nextpmin = paramopts.minAval;
        nextpmax = paramopts.maxAval;
        
    elseif ismember(parameter_names{j}, {'R'})
        
        nextpmin = paramopts.minRval;
        nextpmax = paramopts.maxRval;
        
    elseif ismember(parameter_names{j}, {'S'}) %hack FIXME
        
        nextpmin = paramopts.minSval;
        nextpmax = paramopts.maxSval;
        
    elseif ismember(parameter_names{j}, {'c0'}) %hack FIXME
        
        nextpmin = paramopts.minc0val;
        nextpmax = paramopts.maxc0val;
        
    end
    
    assert(all(nextpmin < nextpmax), 'invalid constraints');
    
    if numel(nextpmin) == 1
        
        nextpmin = repmat(nextpmin, size(pv));
        
    else
        
        assert(ndims(nextpmin) == ndims(pv) && all(size(nextpmin) == size(pv)), 'invalid size');
        
    end
    if numel(nextpmax) == 1
        
        nextpmax = repmat(nextpmax, size(pv));
        
    else
        
        assert(ndims(nextpmax) == ndims(pv) && all(size(nextpmax) == size(pv)), 'invalid size');
        
    end
    
    if any(pv < nextpmin)
        
        fprintf('parameter %s has been raised to its minimum\n', parameter_names{j});
        pv = max(pv, nextpmin);
        
    end
    
    if any(pv > nextpmax)
        
        fprintf('parameter %s has been lowered to its maximum\n', parameter_names{j});
        pv = min(pv, nextpmax);
        
    end
    
    if paramopts.uselog && nextpmin > -inf
        
        pv = log(max(mindiff_init, pv - nextpmin));  % fixme support pmax in uselog case
        
    end
    
    %order values first by element, then by session/fp/etc.    
    p0 = [p0; pv(:)]; %#ok<AGROW>
    pmin = [pmin; nextpmin(:)]; %#ok<AGROW>
    pmax = [pmax; nextpmax(:)]; %#ok<AGROW>
    
    %for each element of the parameter vector, keep a list of segments which depend on it
    for k = 1:ngroups %note that segments with groupindex 0 are not included
        
        for q = 1:nelements
            
            p_seglist_states{1, ii + (k - 1) * nelements + q} = find(groupindex == k); %#ok<AGROW>
            
        end
        
    end
    
    %create indexing arrays to retrieve parameter values from parameter vector for each segment
    pi.(parameter_names{j}) = ii + bsxfun(@plus, (groupindex - 1) * nelements, 1:nelements);
    pi.(parameter_names{j})(groupindex == 0) = NaN;
    
    %build paramter data etc.
    isobs = [isobs; repmat(ismember(parameter_names{j}, paramopts.obsvars), numel(pv), 1)]; %#ok<AGROW>
    
    ii = ii + numel(pv);
end
pi = orderfields(pi);