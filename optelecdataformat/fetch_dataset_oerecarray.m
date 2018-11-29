function oerec = fetch_dataset_oerecarray(oedb, fileinfo, dataset)
if oedb.opt.storedatainmemory
    
    if oedb.issimulation(dataset)
        
        oerec = cat(2, oedb.data{dataset}.oerec);
        
    else
        
        oerec = oedb.data{dataset};
        
    end
    
else
    
    for n = 1:oedb.nneurons(dataset)
        
        oerec(n) = fetch_neurondata(oedb, fileinfo, dataset, n);
        
    end
    
end