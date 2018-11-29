function [oerec_datasets, oesim_datasets, oe_datadir] = oerec_getdatainfo(oe_datadir)
if ~exist('oe_datadir','var')
    oe_datadir = [fileparts(mfilename('fullpath')) filesep '..' filesep 'data' filesep];
end
if ~exist(oe_datadir,'dir')
    oe_datadir = ''; [oerec_datasets, oesim_datasets] = deal({}); return;
end

%use cd to get absolute path without /../ etc.
olddir = cd;
cd(oe_datadir);
oe_datadir = pwd;
cd(olddir);

if oe_datadir(end) ~= filesep
    oe_datadir = [oe_datadir filesep];
end

matfiles_d = dir([oe_datadir '*.mat']);
matfiles = cell(numel(matfiles_d),1);
for u = 1:numel(matfiles_d)
    matfiles{u} = matfiles_d(u).name;
end

oerec_datasets = matfiles(cellfun(@(v) ~isempty(whos('-file',[oe_datadir v],'oerec')), matfiles)); %files containing a variable called oerec
oesim_datasets = matfiles(cellfun(@(v) ~isempty(whos('-file',[oe_datadir v],'oesim')), matfiles)); %files containing a variable called oerec
%should probably run additional checks on each file and spit out warnings for corrupt or improperly formatted files FIXME