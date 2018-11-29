function oedb_cleartmpdir(fileinfo)
if ~isdir(fileinfo.tmpdir)
    
    mkdir(fileinfo.tmpdir);
    
else
    
    delete([fileinfo.tmpdir '*.*']);
    
end