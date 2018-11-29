function r = c2s_predict_wrapper(oerec, modelfile)
inputfilename = [fileparts(mfilename('fullpath')) filesep 'c2s_data_file_tmp.mat'];
outputfilename = [fileparts(mfilename('fullpath')) filesep 'c2s_predict_results_tmp.mat'];
oerec2c2sfile(inputfilename, oerec);

commandstr = sprintf('c2s predict %s %s -v 0', inputfilename, outputfilename);

if exist('modelfile', 'var') && ~isempty(modelfile)
    
    assert(exist(modelfile, 'file') == 2, 'specified modelfile does not exist');
   
    commandstr = sprintf('%s -m %s', commandstr, modelfile);
    
end

assert(system(commandstr) == 0, 'command failed');
r = load(outputfilename, '-mat');
r = r.data;

delete(inputfilename);
delete(outputfilename);