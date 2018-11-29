function oedb = oedb_assign_results(oedb, fileinfo, dataset, neuron, segment, algind, results)
if oedb.opt.storedatainmemory
    
    oedb.results{algind, dataset}{neuron}(segment) = results;
    
else
    
    rfile = rfilename(oedb, fileinfo, algind, dataset, neuron, segment);
    
    if ~isempty(results.spikecounts) || ~isempty(results.spiketimes) || ...
            any(~isnan(results.spiketimes_window)) || ~isempty(fieldnames(results.params)) || ...
            ~isempty(fieldnames(results.opts)) || ~isempty(fieldnames(results.outputvars))
        
        save(rfile, 'results');
        
    elseif exist(rfile,'file')
        
        delete(rfile);
        
    end
    
end