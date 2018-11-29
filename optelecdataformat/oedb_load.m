function [oedb, fileinfo, ui] = oedb_load(filename, fileinfo, prevfigsavedir) %#ok<STOUT>
%[oedb, fileinfo, ui] = oedb_load(filename, fileinfo, prevfigsavedir)
p = fileparts(filename);
if p(end) ~= filesep, p = [p filesep]; end
if ~exist('fileinfo', 'var')
    
    fileinfo = default_oedatabrowser_fileinfostruct;
    
end
if ~exist('prevfigsavedir', 'var')
    
    prevfigsavedir = fileinfo.defaultfigsavedir;
    
end

fileinfo.odbdir = p;
fileinfo.odbfile = filename;

oedb_cleartmpdir(fileinfo);

unzip(fileinfo.odbfile, fileinfo.tmpdir); %unzip works without having to rename to a .zip extension

load([fileinfo.tmpdir 'base.mat'], 'oedb', 'ui'); %overwrites oedb and initializes ui

oedb = update_oedb_version(oedb, fileinfo, oedb_version, min_oedb_version); %#ok<NODEF>

delete([fileinfo.tmpdir 'base.mat']);

if oedb.opt.storedatainmemory %otherwise, necessary data files are in the tmp directory, all ready to go
    
    try %in this try block we don't delete any files, so we can recover and go to streaming mode if we run out of memory
        
        oedb = oedb_loaddatafromdisk(oedb, fileinfo);
        
    catch ex
        
        warning('Failed to load data into memory, will stream from disk instead');
        oedb.opt.storedatainmemory = false;
        oedb = oedb_cleardatafrommemory(oedb);
        
    end
    
end

if oedb.opt.storedatainmemory
    
    oedb_cleartmpdir(fileinfo);
    
end

if ~isdir(oedb.opt.figsavedir)
    
    oedb.opt.figsavedir = prevfigsavedir;
    
end

oedb.modified = false;


function v = min_oedb_version %versions older than this are not supported anymore
v = 0.24;