function [n_isolated_bursts, ndetected] = test_burst_detection(st_true, st_est, tol)
[n_isolated_bursts, ndetected] = deal(0);
if iscell(st_true)
    
    m = numel(st_true);
    for k = 1:m
        
        [ns, nd] = test_burst_detection(st_true{k}, st_est{k}, tol);
        n_isolated_bursts = n_isolated_bursts + ns;
        ndetected = ndetected + nd;
        
    end
    return;
    
end

%define "isolated as > twice tol away from other APs, so there's no ambiguity about which single was matched by an inferred AP
%duration of burst can be up to tol
bursts_starts  = st_true(diff([-inf; st_true(:)]) > 2 * tol & diff([st_true(:); inf]) < tol);
for st = bursts_starts(:)'
    
    last_spike_in_burst = st_true(find(st_true > st & st_true <= st + tol, 1, 'last'));
    if any(st_true > last_spike_in_burst & st_true <= last_spike_in_burst + 2 * tol), continue; end %burst doesn't end within tol
    n_isolated_bursts = n_isolated_bursts + 1;
    if any(st_est >= st - tol & st_est < last_spike_in_burst + tol)
        
        ndetected = ndetected + 1;
        
    end
    
end