function ind = vars_to_indicator(vars_mat,opts)

defaultopts = struct( ...
     'savepath', false ...
);
if ~exist('opts','var')
    opts = struct;
end
opts = mergeOpts(opts,defaultopts);

if numel(vars_mat)>6 && strcmp(vars_mat(1:7),'file://')
    vars_mat = vars_mat(8:end);
end
vars = load(vars_mat);
logSR = log([ unique([vars.Pout.S])' unique([vars.Pout.R])' ] ); % 'unique' is a quick fix to merge multiple segments of one cell. since it is a fitted double value, it is improbable that two will have the same.
                                                                 % TODO: replace with clean solution
if isempty(logSR)
    normal_meancov_logSR = nan(2, 3);  
elseif size(logSR, 1) == 1
    normal_meancov_logSR = [mean(logSR, 1)' nan(2)];    
else
    normal_meancov_logSR = [mean(logSR, 1)' cov(logSR)];
end

ind = struct();
ind.P = struct( ...
     'koff',                    vars.Pout(1).koff ...
    ,'kon',                     vars.Pout(1).kon ...
    ,'A',                       vars.Pout(1).A ...
    ,'c0',                      vars.Pout(1).c0 ...
    ,'tau_ex',                  vars.Pout(1).tau_ex ...
    ,'kd_ex',                   [] ...
    ,'tau_B',                   vars.Pout(1).tau_B ...
    ,'kd_B',                    vars.Pout(1).kd_B ...
    ,'Btot',                    vars.Pout(1).Btot ...
    ,'dbrightness',             vars.Pout(1).dbrightness ...
    ,'normal_meancov_logSR',    normal_meancov_logSR ...
);
ind.V = struct( ...
     'model', '5s2b' ...
    ,'absolute_min_S', 0.4000 ...
); 
ind.uuid = char(java.util.UUID.randomUUID); % TODO: This should probably be saved somewhere in the result file for documentation purposes

if islogical(opts.savepath)
    [f,p] = putfile_fromp(sbm.config_dir, '.mat', 'Save indicator config');
    save(fullfile(p,f),'-mat','-struct','ind');
else
    save(opts.savepath,'-mat'.'-struct','ind');
end

end

