function oedb = oedb_init_outputs(oedb, datasetlist)
if ~exist('datasetlist', 'var')
    
    datasetlist = 1:oedb.ndatasets;
    
end
for algind = 1:numel(oedb.algnames)
    for u = datasetlist(:)'
        
        oedb.results{algind, u} = apdet_resultsstruct(oedb.nsegments{u});
        
    end
end