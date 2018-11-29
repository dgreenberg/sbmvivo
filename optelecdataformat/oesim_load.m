function oesim = oesim_load(filename)
assert(exist(filename,'file') == 2, 'file does not exist');
w = whos('-file',filename);
assert(~isempty(w),'file contains no data');
assert(any(strcmp({w.name},'oesim')),'not a valid simulated optical / electrical recording file');
orig = load(filename,'oesim');
oesim = update_oesim_version(orig.oesim);