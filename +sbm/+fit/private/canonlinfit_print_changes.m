function canonlinfit_print_changes(e0, e, Pfunc, p, p0, parameter_names, pi, R0, R)

%display changes in sum of squared errors and parameters
fprintf('mean ssq error averaged over neurons:    %f -> %f\n',e0,e);

[P0, L0] = Pfunc(p0);
[P,  L] = Pfunc(p);

maxnl = max(cellfun(@numel, parameter_names));
for k = 1:numel(parameter_names)
    
    str = [parameter_names{k} repmat(' ', 1, maxnl - numel(parameter_names{k})) ': '];
    
    if ismember(parameter_names{k}, {'c0' 'kd_B' 'kd_ex' 'A' 'S' 'Btot' 'k50'}) && isfield(P, 'concentration_scale_to_nM')
        
        unitsfac = P(1).concentration_scale_to_nM;
        %fixme check whether concentration_scale_to_nM is consistent?
        
    elseif ismember(parameter_names{k}, {'kon'}) && isfield(P, 'concentration_scale_to_nM')
        
        unitsfac = 1.0 / P(1).concentration_scale_to_nM;  % printed units will be (nM^-1) (s^-1)
        
    else
        
        unitsfac = 1;
        
    end
        
    if ismember(parameter_names{k}, {'dbrightness', 'A', 'kon', 'koff' 'k50'}) && ~isfield(pi, parameter_names{k})
        
        pvals =  unique(cat(1,  P.(parameter_names{k})), 'rows');
        pvals0 = unique(cat(1, P0.(parameter_names{k})), 'rows');
        
    else
        
        piu = unique(pi.(parameter_names{k}));
        piu(isnan(piu)) = [];
        
        if isempty(piu)
        
            continue; %e.g. parameter is not used for this model
            
        else
            
            pvals = L(piu);
            pvals0 = L0(piu);
            
        end
        
    end
    
    for j = 1:numel(pvals0)
        
        str = [str ' ' num2str(reshape(pvals0(j) * unitsfac,1,[]))];
        
    end
    str = [str ' -> '];
    for j = 1:numel(pvals)
        
        str = [str ' ' num2str(reshape( pvals(j) * unitsfac,1,[]))];
        
    end
    fprintf('%s\n',str);
    
end
if exist('R','var') && ~isempty(R)
    
    str = ['R' repmat(' ', 1, maxnl - numel('R')) ': '];
    for j = 1:numel(R0)
        
        str = [str ' ' num2str(R0(j))];
        
    end
    str = [str ' -> '];
    for j = 1:numel(R)
        
        str = [str ' ' num2str(R(j))];
        
    end
    fprintf('%s\n',str);
    
end