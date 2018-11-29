function [dff, tt_states, states] = maxdff_vs_gtot(P, opts, dt, T, gtotlist, ns)
assert(numel(P) == 1, 'P must be a scalar sctruct');
if ~exist('dt','var') || isempty(dt)
    
    dt = 0.005;
    
end
if ~exist('T','var') || isempty(T)
    
    T = 2 + ceil(5 / dt); %fixme bad notation
    
end
if ~exist('gtotlist', 'var')
    
    gtotlist = P.S;
    
end
if ~exist('ns', 'var') || isempty(ns)
    
    ns = 1;
        
end

dff = [];
for gind = 1:numel(gtotlist)
    
    P.S = gtotlist(gind);
    [~, states, tt_states] = sbm.model.spikeresponse(P, opts, dt, T, ns);
    
    if strcmpi(opts.model, '5s2b')
       
        states_eq = states(:, 1);        
        FBGp1 = 1 + P.R * (1 + P.dbrightness(:)' * states_eq(3:6) / P.S);        
        brightness = P.dbrightness(:)' * states(3:6, :) / P.S + FBGp1;
        
    else
        
        error('not yet implemented');
        
    end
    
    dff(:, gind) = 100 * (brightness - brightness(1))' / brightness(1);
    
end