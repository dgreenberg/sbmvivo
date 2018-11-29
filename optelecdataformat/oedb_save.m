function [oedb, fileinfo] = oedb_save(p, f, oedb, fileinfo, ui) %#ok<INUSD>
fileinfo.odbdir = p;
fileinfo.odbfile = [p f];
if p(end) ~= filesep, p = [p filesep]; end

try
    
    datawasinmemory = oedb.opt.storedatainmemory;
    
    if oedb.opt.storedatainmemory %should probably not have two copies of data-full oedb at once, which is what this code effectively does FIXME
        
        oedb = oedb_movedatatodisk(oedb, fileinfo); %sets oedb.opt.storedatainmemory = false;
        
    end
    nfdir = [fileinfo.tmpdir 'newfile' filesep];
    if exist(nfdir, 'dir')
        
        delete([nfdir '*.*']);
        
    else
        
        mkdir(nfdir);
        
    end
    copyfile([fileinfo.tmpdir '*.mat'], nfdir);
    
    save([nfdir 'base.mat'], 'oedb', 'ui'); %note that oedb doesn't currently contain any data, as it's all on disk
    zipfilename = [p f '.zip'];
    zip(zipfilename,'*.mat',nfdir); %zip all mat files in the folder into a single zip file
    delete([nfdir '*.*']);
    rmdir(nfdir);
    movefile(zipfilename, [p f]); %the zip function adds a .zip extension even if another is specified
    
    if datawasinmemory
        
        oedb = oedb_movedatatomemory(oedb, fileinfo);
        
    end
    oedb.modified = false;
    
catch ex
    
    uiwait(errordlg(ex.message, 'Failed to save file', 'modal')); return; %note that in case of a crash, we haven't overwritten oedb with a setappdata call
    %fixme figure out whether the program state is now invalid, why we
    %failed etc., modify fileinfo
    
end