function [P_eachseg, L] = pvec2Pstruct(P0, parameter_names, pi, p, pmin, opts)
L = lindomain(p, pmin, opts);

P_eachseg = P0;

for j = 1:numel(parameter_names)
    
    nexti = pi.(parameter_names{j}); %nexti is a vector of indices into p/L for each segment, for this parameter
    segments_used = find(~any( isnan(nexti), 2));
    for k = reshape(segments_used, 1, []) %skip segments that were not used due to lack of baseline or spikes
        
        P_eachseg(k).(parameter_names{j}) = reshape(L(nexti(k,:)), 1, []);
        
    end
    
end

P_eachseg = update_dependent_params_autoR(P_eachseg, opts, segments_used);


function L = lindomain(p, pmin, paramopts)
L = p;
if paramopts.uselog
    
    L(pmin > -inf) = exp(p(pmin > -inf)) + pmin(pmin > -inf);
    
end
