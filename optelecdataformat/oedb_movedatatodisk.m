function oedb = oedb_movedatatodisk(oedb, fileinfo)
assert(oedb.opt.storedatainmemory, 'attempted to move data to disk, but data is already there');
oedb_cleartmpdir(fileinfo);
oedb_savedatatodisk(oedb, fileinfo);
oedb = oedb_cleardatafrommemory(oedb);
oedb.opt.storedatainmemory = false;