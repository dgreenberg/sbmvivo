function results = fetch_results(oedb, fileinfo, dataset, neuron, segment, algind)
%results = fetch_results(oedb, fileinfo, dataset, neuron, segment, algind)
nseg = 1;

if oedb.opt.storedatainmemory
    
    results = oedb.results{algind, dataset}{neuron}(segment);
    
else
    
    results = apdet_resultsstruct(nseg); %default empty struct to return if the data is not present
    results = results{1};
    
    rfile = rfilename(oedb, fileinfo, algind, dataset, neuron, segment);
    
    if exist(rfile, 'file')
        
        w = warning('query','MATLAB:load:variableNotFound');
        warning('off','MATLAB:load:variableNotFound'); %file can be missing a variable due to old file version etc., no need to spit out warning in that case
        load([fileinfo.tmpdir 'r' num2str(algind) '_' oedb.segcode{dataset}{neuron}{segment} '.mat'], 'results');
        warning(w.state, 'MATLAB:load:variableNotFound');
        
    end
    
end