function oedb = oedb_cleardatafrommemory(oedb)
oedb.data = cell(1, oedb.ndatasets);
for u = 1:oedb.ndatasets
    for v = 1:oedb.nneurons(u)
        for algind = 1:numel(oedb.algnames)
            
            oedb.results{algind, u} = cell(1, oedb.nneurons(u));
            
        end
    end
end