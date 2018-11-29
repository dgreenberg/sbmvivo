function oedb = oedb_init_stats(oedb, datasetlist)
if ~exist('datasetlist', 'var')
    
    datasetlist = 1:oedb.ndatasets;
    
end
nalgs = numel(oedb.algnames);

for u = datasetlist(:)'
    
    oedb.stats.bydataset(:, u) = repmat(empty_statstruct, nalgs, 1);
    
    oedb.stats.byneuron{u} = repmat(empty_statstruct, nalgs, oedb.nneurons(u));
    
    oedb.stats.bysegment{u} = cell(1, oedb.nneurons(u));
    
    for n = 1:oedb.nneurons(u)
        
        oedb.stats.bysegment{u}{n} = repmat(empty_statstruct, nalgs, oedb.nsegments{u}(n));
        
    end
    
end