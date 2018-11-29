function [npairs, sa, sb, amissed, bmissed, saok_any, sbok_any] = test_spike_recon(sta, stb, tol, wa, wb)
%[npairs, na, nb, amissed, bmissed, saok_any, sbok_any] = test_spike_recon(sta, stb, tol, wa, wb)
if iscell(sta) %multiple spike train pairs
    
    m = numel(sta);
    if ~exist('wa','var'), wa = cell(1, m); end
    if ~exist('wb','var'), wb = cell(1, m); end      
    
    [npairs, sa, sb, amissed, bmissed, saok_any, sbok_any] = deal(nan(1, m));
    
    for k = 1:m
        
        [npairs(k), sa(k), sb(k), amissed(k), bmissed(k), saok_any(k), sbok_any(k)] = test_spike_recon(sta{k}, stb{k}, tol, wa{k}, wb{k});        
        
    end
    return;
    
end
%single spike train:
na = numel(wa); nb = numel(wb);
[npairs, saok_any, sbok_any, amissed, bmissed] = deal(0);
if ~exist('wa','var') || isempty(wa), wa = ones(size(sta)); end
if ~exist('wb','var') || isempty(wb), wb = ones(size(stb)); end

sa = sum(wa); sb = sum(wb);

if ~na || ~nb, return; end
assert(size(sta,2) == 1 && size(stb,2) == 1,'sta and stb must be column vectors');
assert(all(diff(sta) >= 0) && all(diff(stb) >= 0), 'spike time lists must be pre-sorted');
assert(all(wa > 0) && all(wb > 0), 'weights must be positive');

tdiff = sta * ones(1, nb) - ones(na, 1) * stb';
nears = sparse(abs(tdiff) <= tol);
saok_any = sum(wa(any(nears,2)));
sbok_any = sum(wb(any(nears,1)));
if ~any(any(nears))
    
    return;
    
end
todoa = any(nears,2);
todob = any(nears,1);
npairs = 0;
while(any(todoa) && any(todob))
    
    todoa = any(nears(:,todob),2) & todoa; 
    todob = any(nears(todoa,:),1) & todob;
    ia = find(todoa, 1);
    ib = find(todob, 1);
    
    %let the weights (default 1) cancel out between the two spike trains:
    dw = min(wa(ia), wb(ib));
    wa(ia) = wa(ia) - dw;
    if wa(ia) == 0
        
        todoa(ia) = 0; %this spike is totally canceled
        
    end
    wb(ib) = wb(ib) - dw;
    if wb(ib) == 0
        
        todob(ib) = 0; %this spike is totally canceled
        
    end
    
    npairs = npairs + dw; %keep track of total weight of matched spikes
    
end
amissed = sa - npairs;
bmissed = sb - npairs;