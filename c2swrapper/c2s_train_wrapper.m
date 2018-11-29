function c2s_train_wrapper(oerec, modelfilename)
inputfilename = [fileparts(mfilename('fullpath')) filesep 'c2s_data_file_tmp.pck'];
oerec2c2sfile(inputfilename, oerec);

commandstr = sprintf('c2s train %s %s -v 0', inputfilename, modelfilename);
assert(system(commandstr) == 0, 'command failed');