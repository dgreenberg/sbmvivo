function oedb = oedb_movedatatomemory(oedb, fileinfo)
assert(~oedb.opt.storedatainmemory, 'attempted to move data to memory, but data is already there');
oedb = oedb_loaddatafromdisk(oedb, fileinfo);
oedb_cleartmpdir(fileinfo);
oedb.opt.storedatainmemory = true;