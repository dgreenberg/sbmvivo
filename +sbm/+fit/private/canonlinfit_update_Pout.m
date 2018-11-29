function [Pout, Rind] = canonlinfit_update_Pout(Pout, focalplaneindex, neuronindex, focalplaneindex_orig, neuronindex_orig, opts, R)
%assign R values
nseg = numel(Pout);
if ismember('R', opts.focalplanespecific)
    
    Rind = focalplaneindex;
    Rind_orig = focalplaneindex_orig;
    
elseif ismember('R', opts.neuronspecific)
    
    Rind = neuronindex;
    Rind_orig = neuronindex_orig;
    
else
    
    [Rind, Rind_orig] = deal(ones(nseg, 1));
    
end
for k = 1:nseg
    if Rind(k) > 0 %excludes NaN
        
        Pout(k).R = R(Rind(k));
        
    else  % we didn't use this segment because for example we couldn't find baseline or the segment was too short
        
        ii = Rind_orig == Rind_orig(k) & Rind > 0;
        if any(ii)
            
            Pout(k).R = mean(R(Rind(ii)));
            
        else
            
            Pout(k).R = nan;
            
        end
    end
end
for u = 1:numel(opts.neuronspecific)
    for k = 1:nseg
        
        if any(isnan(Pout(k).(opts.neuronspecific{u})))
            
            ii = neuronindex(neuronindex_orig == neuronindex_orig(k) & neuronindex > 0);
            vals = cat(1, Pout(ii).(opts.neuronspecific{u}));
            ok = ~isnan(vals);
            vals(~ok) = 0;
            Pout(k).(opts.neuronspecific{u}) = sum(vals, 1) ./ sum(ok, 1);
            
        end
    end
end
for u = 1:numel(opts.focalplanespecific)
    for k = 1:nseg
        
        if any(isnan(Pout(k).(opts.focalplanespecific{u})))
            
            ii = focalplaneindex(focalplaneindex_orig == focalplaneindex_orig(k) & focalplaneindex > 0);
            vals = cat(1, Pout(ii).(opts.focalplanespecific{u}));
            ok = ~isnan(vals);
            vals(~ok) = 0;
            Pout(k).(opts.focalplanespecific{u}) = sum(vals, 1) ./ sum(ok, 1);
            
        end
    end
end
if ~opts.usek50
    
    for k = 1:nseg
        
        [Pout(k).k50, Pout(k).tau_mm] = koffkon2k50tau(Pout(k).koff, Pout(k).kon);
        Pout(k).k50_first = Pout(k).k50(1);
        Pout(k).dk50 = diff(Pout(k).k50(1));
        
    end
    
end