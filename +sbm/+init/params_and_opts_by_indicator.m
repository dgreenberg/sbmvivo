function [P,V] = params_and_opts_by_indicator(indicatorstring)
% Function imports mat files with fitting parameters as two structs of the 
% following shape
%   V = struct( ...
%       'model, '5s2b' ...
%      ,'absolute_min_S, 0.4 ... % uM
%       );
%
%   P = struct( ...
%        'koff',  [0.1001686918 204.7802878 11.77442605 5.831340864] ...  %fixme: Document values
%       ,'kon',    [2.516340108  16.93083934 1.087582073 1068.893443] ...
%       ,'A',      20.1697 ...
%       ,'c0',     0.0500 ...
%       ,'tau_ex', 0.007489 ...
%       ,'kd_ex', [] ...
%       ,'tau_B',  [0.001       17.11292069] ...
%       ,'kd_B',   [3.391135289 0.5946177344] ...
%       ,'Btot',   [63.77391739 118.9676718] ...
%       ,'dbrightness',  [0 0 0 80] ...
%       ,'normal_meancov_logSR', [1.9996    6.1610   -3.0240; 0.9347   -3.0240    5.9039] ...  % joint normal prior over logs of background fluorescence and indicator concentrations
%       );
%
% normal_meancov_logSR is calculated as
%   logSR = log([S_eachneuron, R_eachneuron]); % _eachneuron has shape n x 1
%     or
%   logSR = log([ unique([Pout.S])' unique([Pout.R])' ] )
%   normal_meancov_logSR = [mean(logSR, 1)' cov(logSR)];
% All other quantities in P are direct outputs of oedb_nonlinfit.m

indicatorstring = normalize_indicator_string(indicatorstring);

files = { ...
    fullfile(sbm.config_dir, 'indicators', 'ind_%s.mat'), ...
    fullfile(fileparts(mfilename('fullpath')), '..', 'indicators', 'ind_%s.mat') ...
};
% fixme: There should be third method where people provide a user defined
% directory to allow group-internal repositories or network drives.

files = cellfun(@(s) sprintf(s,indicatorstring),files,'UniformOutput',false);

for i = 1:numel(files)
    
    if exist(files{i},'file')
    
        fprintf('Using file %s for indicator parameters.\n', files{i});
        load(files{i},'-mat','P','V');
        break;
        
    end
    
end

if ~exist('P','var') || ~exist('V','var')
    %fixme: we should suppress this error if all needed params were supplied
    %by the user, i.e., check for all required fields
    
    error('sbm:indicatornotrecognized','unrecognized indicator');
        
end

end

function indicatorstring = normalize_indicator_string(indicatorstring)
% TODO: If necessary (as indicated by old version of main function)
% this should probably load csv-files from repository and user
% directory for string translation. 
indicatorstring = lower(indicatorstring);
end