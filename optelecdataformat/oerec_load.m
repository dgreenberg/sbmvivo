function oerec = oerec_load(filename)
assert(exist(filename,'file') == 2, 'file does not exist');
w = whos('-file',filename);
assert(~isempty(w),'file contains no data');
assert(any(strcmp({w.name},'oerec')),'not a valid optical / electrical recording file');
orig = load(filename,'oerec');
oerec = update_oerec_version(orig.oerec);