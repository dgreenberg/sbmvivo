function [n_isolated_singles, ndetected] = test_singleAP_detection(st_true, st_est, tol)
[n_isolated_singles, ndetected] = deal(0);
if iscell(st_true)
    
    m = numel(st_true);
    for k = 1:m
        
        [ns, nd] = test_singleAP_detection(st_true{k}, st_est{k}, tol);
        n_isolated_singles = n_isolated_singles + ns;
        ndetected = ndetected + nd;
        
    end
    return;
    
end

%define "isolated as > twice tol away from other APs, so there's no ambiguity about which single was matched by an inferred AP
isolated_singles  = st_true(diff([-inf; st_true(:)]) > 2 * tol & diff([st_true(:); inf]) > 2 * tol);
n_isolated_singles = numel(isolated_singles);
for st = isolated_singles(:)'
    
    if any(st_est >= st - tol & st_est < st + tol)
        
        ndetected = ndetected + 1;
        
    end
    
end