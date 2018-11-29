function print_param_changes(P, V, params_to_print)
ndigits = 6;
fprintf('\n');
if ~exist('params_to_print','var') || isempty(params_to_print)
    
    params_to_print = {'gain' 'zeta' 'fr' 'S' 'R'}; %'sigma_r'
    
end

show_percent_params = {'gain' 'zeta' 'fr' 'S' 'R'}; %'sigma_r'
[fn_disp, maxnamelength] = param_display_names(params_to_print, V);

for pi = 1:numel(params_to_print)
    
    fn = params_to_print{pi};
    pv1 = cat(1, P(end-1,:).(fn))';
    pv2 = cat(1, P(end,:).(fn))';
    
    if size(pv1, 2) == 1
        
        pvboth = [pv1 pv2];
        pvbothunique = unique(pvboth, 'rows', 'stable');
        pv1 = pvbothunique(:, 1);
        pv2 = pvbothunique(:, 2);
        
    else
        
        if issame_ndigits(pv1, ndigits)
            
            pv1 = pv1(:,1);
            
        end
        
        if issame_ndigits(pv2, ndigits)
            
            pv2 = pv2(:,1);
            
        end
        
    end
    
    if issame_ndigits([pv1 pv2], ndigits)
        
        changemark = '===';
        
    else
        
        changemark = '-->';
        
    end
    
    fprintf('%s: ' ,[fn_disp{pi} repmat(' ',1,max(0, maxnamelength - numel(fn_disp{pi})))]);
    
    print_numeric_val(pv1, ndigits, 12);
    
    fprintf('   %s  ',changemark);
    
    print_numeric_val(pv2, ndigits, 12);
    
    if ismember(fn, show_percent_params) && numel(pv1) == 1 && numel(pv2) == 1 && ~issame_ndigits([pv1 pv2], ndigits)
        
        percent_change = 100 * (pv2 - pv1) / pv1;
        valstring = num2str(percent_change, 3);
        if percent_change >= 0, valstring = ['+' valstring]; end %#ok<AGROW>
        fprintf('  %s%%', valstring);
        
    end
    
    fprintf('\n');
    
end


function s = issame_ndigits(pv, ndigits) %see if all elements of a vector will be displayed identically on scren given ndigits decimal places
s = true;
for u = 2:size(pv,2)
    for v = 1:size(pv,1)
        if ~strcmp(num2str(pv(v, 1), ndigits), num2str(pv(v, u), ndigits))
            s = false; return;
        end
    end
end


function [fn_disp, maxnamelength] = param_display_names(params_to_print, V)
fn_disp = cell(size(params_to_print));
for pi = 1:numel(params_to_print)
    
    fn = params_to_print{pi};
    fn_disp{pi} = fn;
    
end
maxnamelength = max(cellfun(@numel, fn_disp));