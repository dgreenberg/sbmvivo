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
        
    case 'gcamp'
        
        if ~strcmpi(gcampsuffix, '6s')
            
            error('unrecognized gcamp variant: %s', gcampsuffix);
            
        end
        
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
        
    
    otherwise
        
        %fixme: we should suppress this error if all needed params were supplied
        %by the user.
        error('sbm:indicatornotrecognized','unrecognized indicator');
        
end