function r = c2s_loo_wrapper_1tn(oerec)
%[ns, r] = c2s_loo_wrapper_1tn(oerec)
%run c2s with 1 training neuron (always the next neuron, looping around)
nneurons = numel(oerec);
ns = cell(1, nneurons);
inputfilename = [fileparts(mfilename('fullpath')) filesep 'c2s_data_file_tmp.mat'];
outputfilename = [fileparts(mfilename('fullpath')) filesep 'c2s_loo_results_tmp.mat'];
r = struct('data', {{}});
for n = 1:nneurons
    
    training_neuron = mod(n, nneurons) + 1;
    oerec_sub = oerec([n training_neuron]);
    
    oerec2c2sfile(inputfilename, oerec_sub);    
    commandstr = sprintf('c2s leave-one-out %s %s', inputfilename, outputfilename);
    assert(system(commandstr) == 0, 'command failed');    
    
    rnext = load(outputfilename, '-mat');
    rnext.data = rnext.data(1:numel(oerec_sub(1).data));
    r.data = [r.data rnext.data];
    
    delete(inputfilename);
    delete(outputfilename);
    
end