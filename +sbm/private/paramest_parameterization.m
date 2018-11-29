function [Pfunc, pvec0, pnames, index] = paramest_parameterization(params0, opts)
nseg = numel(params0);

[~, ~, settingsindex_renumbered] = unique(opts.settingsindex);
n_settings = max(settingsindex_renumbered);
[~, ~, sessionindex_renumbered] = unique(opts.sessionindex);
n_sessions = max(sessionindex_renumbered);
[~, ~, focalplaneindex_renumbered] = unique(opts.focalplaneindex);
n_focalplanes = max(focalplaneindex_renumbered);

pnames = {'S', 'R'};
if ~opts.profile_firing_rate
    %pnames{1, end + 1} = 'fr';
end
if opts.shotnoise
    pnames{1, end + 1} = 'gain';
end
if (opts.darknoise && ~opts.darknoisemeasured) && ~opts.profile_zeta
    %pnames{1, end + 1} = 'zeta';
end
if ~opts.darknoise && ~opts.shotnoise %FIXME are these supposed to indicate model class, what we know from calibration, which parameters we want to fix, or what?
    warning('neither shotnoise nor dark noise are enabled, data-fitting is likely to fail!');
end

if ~opts.darkoffsetcorrected
    pnames{1, end + 1} = 'fdc';
end

%FIXME what kind of parameter should sigma_r be? day-specific perhaps?
neuronspecific = {'S', 'fr', 'R'}; %S,R should probably be day-specific too
settingsspecific = {'gain', 'zeta', 'fdc'};
sessionspecific = {};
focalplanespecific = {};
positive_params = {'S', 'gain', 'zeta', 'sigma_r', 'R' 'fr'};
unitinterval_params = {};
index = cell(1, numel(pnames));
nfreeparams = 0;
for j = 1:numel(pnames)
    
    P0vals_were_averaged = false;
    if ismember(pnames{j}, neuronspecific)
        
        P0vals = cat(1, params0.(pnames{j}));
        if any(any(diff(P0vals, 1, 1) ~= 0))
            
            P0vals_were_averaged = true;
            
        end
        assert(size(P0vals, 2) == 1, 'real-valued params only');
        pvec0(nfreeparams + 1, 1:size(P0vals, 2)) = Pfield2vecelements(mean(P0vals, 1), pnames{j}, positive_params, unitinterval_params); %#ok<AGROW>
        index{j} = ones(nseg, 1) * (nfreeparams + 1);
        nfreeparams = nfreeparams + 1;
        
    elseif ismember(pnames{j}, focalplanespecific)
        
        for k = 1:n_focalplanes
            
            P0vals = cat(1, params0(focalplaneindex_renumbered == k).(pnames{j}));
            if any(any(diff(P0vals, 1, 1) ~= 0))
                
                P0vals_were_averaged = true;
                
            end
            assert(size(P0vals, 2) == 1, 'real-valued params only');
            pvec0(nfreeparams + k, 1:size(P0vals, 2)) = Pfield2vecelements(mean(P0vals, 1), pnames{j}, positive_params, unitinterval_params);
            
        end
        index{j} = nfreeparams + focalplaneindex_renumbered;
        nfreeparams = nfreeparams + n_focalplanes;
        
    elseif ismember(pnames{j}, settingsspecific)
        
        for k = 1:n_settings
            
            P0vals = cat(1, params0(settingsindex_renumbered == k).(pnames{j}));
            if any(any(diff(P0vals, 1, 1) ~= 0))
                
                P0vals_were_averaged = true;
                
            end
            assert(size(P0vals, 2) == 1, 'real-valued params only');
            pvec0(nfreeparams + k, 1:size(P0vals, 2)) = Pfield2vecelements(mean(P0vals, 1), pnames{j}, positive_params, unitinterval_params);
            
        end
        index{j} = nfreeparams + settingsindex_renumbered;
        nfreeparams = nfreeparams + n_settings;
        
    elseif ismember(pnames{j}, sessionspecific)
        
        for v = 1:n_sessions
            
            P0vals = cat(1, params0(sessionindex_renumbered == v).(pnames{j}));
            if any(any(diff(P0vals, 1, 1) ~= 0))
                
                P0vals_were_averaged = true;
                
            end
            assert(size(P0vals, 2) == 1, 'real-valued params only');
            pvec0(nfreeparams + v, 1:size(P0vals, 2)) = Pfield2vecelements(mean(P0vals, 1), pnames{j}, positive_params, unitinterval_params);
            
        end
        index{j} = nfreeparams + sessionindex_renumbered;
        nfreeparams = nfreeparams + n_sessions;
        
    else
        
        error('parameter type not specified for %s', pnames{j});
        
    end
    if P0vals_were_averaged
        
        warning('multiple values for %s were averaged to create initial parameter vector', pnames{j});
        
    end
    
end
Pfunc = @(pvec) segmentparams(params0, pnames, index, pvec, positive_params, unitinterval_params);


function vv = Pfield2vecelements(fv, pname, positive_params, unitinterval_params)
if ismember(pname, positive_params)
    
    vv = log(fv);
    
elseif ismember(pname, unitinterval_params)
    
    vv = log(fv ./ (1 - fv));
    
else
    
    vv = fv;
    
end


function fv = vecelements2Pfield(vv, pname, positive_params, unitinterval_params)
if ismember(pname, positive_params)
    
    fv = exp(vv);
    
elseif ismember(pname, unitinterval_params)
    
    fv = exp(vv) ./ (1 + exp(vv));
    
else
    
    fv = vv;
    
end


function P = segmentparams(P0, pnames, index, pvec, positive_params, unitinterval_params)
P = P0;
for j = 1:numel(pnames)
    
    for k = 1:numel(P0)
        
        if ~index{j}(k)
            
            continue;
            
        end %this parameter is not a free variable for this particular segment (e.g. dark offset was measured for some but not all data)
        P(k).(pnames{j}) = vecelements2Pfield(reshape(pvec(index{j}(k)), 1, []), pnames{j}, positive_params, unitinterval_params);
        
    end
    
end