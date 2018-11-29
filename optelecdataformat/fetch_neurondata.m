function oerec = fetch_neurondata(oedb, fileinfo, dataset, neuron)
%oerec = fetch_neurondata(oedb, fileinfo, dataset, neuron)
if oedb.issimulation(dataset)
    
    oesim = fetch_simstruct(oedb, fileinfo, dataset, neuron);
    oerec = oesim.oerec;
    
elseif oedb.opt.storedatainmemory
    
    oerec = oedb.data{dataset}(neuron);
    
else
    
    oerec = oerec_load([fileinfo.tmpdir 'n' oedb.neuroncode{dataset}{neuron} '.mat']); %also updates version if necessary
    
end