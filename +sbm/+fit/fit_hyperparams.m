function P = fit_hyperparams(P, neuronind, firingrates_eachneuron)
nlist = unique(neuronind(:)');
nneurons = numel(nlist);
[R_eachneuron, S_eachneuron] = deal(nan(nneurons, 1));
for k = 1:nneurons
    
    ii = find(neuronind == nlist(k), 1);
    if isempty(ii) || isempty(P(ii).R), continue; end  % e.g. no baseline available for this neuron
    R_eachneuron(k) = P(ii).R;
    S_eachneuron(k) = P(ii).S;
    
end
logSR = log([S_eachneuron, R_eachneuron]);
logSR(any(isnan(logSR), 2), :) = [];
if any(R_eachneuron == 0) || any(S_eachneuron == 0)
    
    warning('S or R is zero, cannot fit hyperparamters');
    
end
if isempty(logSR)
    
    normal_meancov_logSR = nan(2, 3);
    
elseif size(logSR, 1) == 1
    
    normal_meancov_logSR = [mean(logSR, 1)' nan(2)];
    
else
    
    normal_meancov_logSR = [mean(logSR, 1)' cov(logSR)];
    
end
for j = 1:numel(P)
    
    P(j).normal_meancov_logSR  = normal_meancov_logSR;
    
end
if nargin > 2
    
    fr = firingrates_eachneuron(~isnan(firingrates_eachneuron));
    kt = gamma_fit_kt(fr);
    for j = 1:numel(P)
        
        P(j).gamma_kt_FR = kt;
        
    end
    
end