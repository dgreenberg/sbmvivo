function P_eachseg = update_dependent_params_autoR(P_eachseg, opts, segments_used)
if opts.usek50
    
    [k50_prev, tau_mm_prev] = deal(nan(1, opts.nbindingsteps));
    for k = reshape(segments_used, 1, [])
        
        k50 = cumsum([P_eachseg(k).k50_first P_eachseg(k).dk50]);
        P_eachseg(k).k50 = k50;
        if all(k50 == k50_prev) && all(P_eachseg(k).tau_mm == tau_mm_prev)
            
            [P_eachseg(k).koff, P_eachseg(k).kon] = deal(koff_prev, kon_prev); %use cached copy for speed
            
        else
            
            [P_eachseg(k).koff, P_eachseg(k).kon] = k50tau2koffkon(k50, P_eachseg(k).tau_mm);
            k50_prev = k50;
            tau_mm_prev = P_eachseg(k).tau_mm;
            [koff_prev, kon_prev] = deal(P_eachseg(k).koff, P_eachseg(k).kon);
            
        end
        
    end
    
end
for k = reshape(segments_used, 1, [])
    
    frac_of_max = cumprod(P_eachseg(k).brightness_ratio(end:-1:1));
    frac_of_max = frac_of_max(end:-1:1);
    
    P_eachseg(k).dbrightness = P_eachseg(k).fratio * [frac_of_max 1];
    
end