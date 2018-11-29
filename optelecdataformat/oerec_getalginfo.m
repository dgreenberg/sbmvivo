function [oerec_algfunctions, oerec_alginfo, oerec_algdir, oerec_algfiles, oerec_trainfunctions] = oerec_getalginfo(warnmissing)
%[oerec_algfunctions, oerec_alginfo, oerec_algdir, oerec_algfiles, oerec_trainfunctions] = oerec_getalginfo
%
%algs take a three inputs:
%1) an oerec struct array over neurons. see empty_oerec.m
%2) a cell array (over neurons) of cell arrays (over segements) of structs of parameters
%3) a struct of options
%
%algs must have the following outputs:
%1) results: a cell array over neurons, with each cell a struct array of data segments
%2) an algorithm info struct
if ~exist('warnmissing', 'var') || isempty(warnmissing), warnmissing = false; end

oerec_algdir = [fileparts(mfilename('fullpath')) filesep '..' filesep 'algorithms' filesep];

if ~exist(oerec_algdir, 'dir')
    
    oerec_algdir = '';
    [oerec_algfiles, oerec_algfunctions, oerec_trainfunctions] = deal({});
    oerec_alginfo = empty_apdetalginfo;
    oerec_alginfo(:) = [];
    warning('oerec:noalgdir','could not locate algorithms directory');
    return;
    
end

olddir = cd;
cd(oerec_algdir); %in case the function appears elsewhere on the matlab path
oerec_algfiles_d = dir([oerec_algdir 'apdet_*.m']);
oerec_algfiles = reshape({oerec_algfiles_d.name}, [], 1);
oerec_algfiles = oerec_algfiles(cellfun(@(v) nargin([oerec_algdir v]) == 3, oerec_algfiles) & cellfun(@(v) nargout([oerec_algdir v]) == 2, oerec_algfiles));

oerec_alginfo = repmat(empty_apdetalginfo, 1, numel(oerec_algfiles));
[oerec_algfunctions, oerec_trainfunctions] = deal(cell(1, numel(oerec_algfiles)));
badfiles = false(1, numel(oerec_alginfo));

%empty oerec struct to pass when getting info for each algorithm
oerec = empty_oerec;
oerec(:) = [];
[params, opts] = deal(struct());

for u = 1:numel(oerec_alginfo)
    
    [~, funcname, ~] = fileparts(oerec_algfiles{u});
    oerec_algfunctions{u} = str2func(funcname);
    
    try
        
        [~, oerec_alginfo(u)] = oerec_algfunctions{u}(oerec, params, opts);
        
    catch ex
        
        badfiles(u) = true;
        if warnmissing
            
            warning('oerec:badalginfocall','failed to get algorithm info for %s: %s', oerec_algfiles{u}, ex.message);
            
        end
        
    end
    
    trainfile = ['train_' funcname(7:end)];
    if exist(trainfile, 'file') == 2
        
        oerec_trainfunctions{u} = str2func(trainfile);
        
    end
    
end
oerec_algfunctions(badfiles) = []; oerec_alginfo(badfiles) = []; oerec_algfiles(badfiles) = []; oerec_trainfunctions(badfiles) = [];
cd(olddir);

%put sbm first:
sbmind = find(strcmpi({oerec_alginfo.shortname}, 'sbm'));
if numel(sbmind) == 1
    
    neworder = [sbmind; setdiff((1:numel(oerec_alginfo))', sbmind)];
    
    oerec_algfunctions   = oerec_algfunctions(neworder);
    oerec_alginfo        = oerec_alginfo(neworder);
    oerec_algfiles       = oerec_algfiles(neworder);
    oerec_trainfunctions = oerec_trainfunctions(neworder);
    
end