function oedb = oedb_clear_results(oedb, fileinfo, handles, algind, datasetind, neuronind, segmentind)
assert(numel(datasetind) == 1, 'invalid input');
results = apdet_resultsstruct(1); %initialize results struct for one datasegment
results = results{1};

if ~exist('neuronind','var') || isempty(neuronind) %entire dataset
    
    assert(~exist('segmentind','var') || isempty(segmentind), 'segment can''t be specified unless neuron is specified');
    
    oedb.stats.bydataset(algind, datasetind) = empty_statstruct;
    
    oedb.stats.byneuron{datasetind}(algind, :) = repmat(empty_statstruct, 1, oedb.nneurons(datasetind));
    
    for ni = 1:oedb.nneurons(datasetind)
        
        oedb.stats.bysegment{datasetind}{ni}(algind, :) = repmat(empty_statstruct, 1, oedb.nsegments{datasetind}(ni));
        
        for si = 1:oedb.nsegments{datasetind}(ni)
            
            oedb = oedb_assign_results(oedb, fileinfo, datasetind, ni, si, algind, results);
            
        end
        
    end
    
else
    
    oedb.stats.byneuron{datasetind}(algind, neuronind) = repmat(empty_statstruct, 1, numel(neuronind));
    
    for ni = neuronind
        
        oedb.stats.bysegment{datasetind}{ni}(algind, :) = repmat(empty_statstruct, 1, oedb.nsegments{datasetind}(ni));
        
        if ~exist('segmentind','var') || isempty(segmentind)
            
            segmentind = 1:oedb.nsegments{datasetind}(ni);
            
        end
        
        for si = segmentind
            
            oedb = oedb_assign_results(oedb, fileinfo, datasetind, ni, si, algind, results);
            
        end
        
        oedb = oedb_calc_neuron_stats(oedb, fileinfo, algind, datasetind, ni);
        
    end
    
end