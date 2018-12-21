function [P,V] = params_and_opts_by_indicator(indicatorstring)
gcloc = strfind(lower(indicatorstring), lower('GCaMP'));
if ~isempty(gcloc)
    
    if numel(indicatorstring) >= gcloc + 5
        gcampsuffix = indicatorstring(gcloc + 5:end);
        gcampsuffix = gcampsuffix(find(gcampsuffix ~= ' ',1):end);
    end
    indicatorstring = 'gcamp';
    
end
switch lower(indicatorstring) %assign indicator-dependent parameters
    case 'ogb1-am'
        
        V.onrate = false;
        V.saturation = false;
        V.poly = false;
        P.mu_b_init = 50;
        P.c0      = 25;
        P.A       = 250;
        P.k_d     = 300;
        P.tau_c   = 0.5;
        P.n       = 1; %hill exponent
        V.expected_Frat_per_ap = 0.1; %expected frac. increase in F from 1 AP
        V.bound_to_free_fluorescence_ratio = 14; %"Long-Wavelength Calcium Indicators" by Molecular Probes http://tools.invitrogen.com/content/sfs/manuals/mp03010.pdf
        
    case 'gcamp'
        
        V.model = '5s2b';
        if strcmpi(V.model, '5s2b')
            
            %micromolar units
            
            P.koff = [0.1001686918 204.7802878 11.77442605 5.831340864];
            P.kon =  [2.516340108  16.93083934 1.087582073 1068.893443];
            
            P.A = 20.1697;
            P.c0 = 0.0500;
            P.tau_ex =  0.007489;
            P.kd_ex  = [];
            
            P.tau_B = [0.001       17.11292069];
            P.kd_B =  [3.391135289 0.5946177344];
            P.Btot =  [63.77391739 118.9676718];
            
            P.dbrightness =  [0 0 0 80];
            
            P.normal_meancov_logSR = [1.9996    6.1610   -3.0240; 0.9347   -3.0240    5.9039];  % joint normal prior over logs of background fluorescence and indicator concentrations
            
            V.absolute_min_S = 0.4; % uM
            
        else
            
            error('unrecognized model: %s', V.model);
            
        end
        
        if gcampsuffix(1) == '3'
            
            P.n = 2.54;
            P.k_d = 345 ^ P.n;
            P.tau_s = 1 / 2.57;
            
        elseif gcampsuffix(1) == '5' %for 5G, not sure about other variants
            
            P.tau_s = 1 / 2.52;
            
        elseif gcampsuffix(1) == '6'
            
            if strcmpi(gcampsuffix(2),'m')
                
                V.expected_Frat_per_ap = 0.13;
                
            elseif strcmpi(gcampsuffix(2),'f')
                
                V.expected_Frat_per_ap = 0.19;
                
            elseif strcmpi(gcampsuffix(2),'s') || strcmpi(gcampsuffix,'641')
                
                V.expected_Frat_per_ap = 0.14; %0.23;
                
            else
                
                error('unrecognized GCaMP6 variant');
                
            end
            
        else
            
            error('sbm:gcampvariantnotrecognized','unrecognized GCaMP variant');
            
        end
        
    otherwise
        
        %fixme: we should suppress this warning if all params were supplied
        %by the user.
        error('sbm:indicatornotrecognized','unrecognized indicator');
        V.onrate = false; %FIXME don''t copy-paste these values, it'll lead to a mistake later
        
end