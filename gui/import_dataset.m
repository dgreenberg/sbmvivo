function oedb = import_dataset(oedb, fileinfo, datafile)
assert(exist(datafile, 'file') == 2, 'file not found');

u = oedb.ndatasets + 1;
oedb.sourcedatafiles{u} = datafile;

nneuronstotal = sum(oedb.nneurons);
nsegtotal = sum(cat(2, oedb.nsegments{:}));

[~, ff, ee] = fileparts(datafile);
if numel(ff) > 0 && strcmpi(ee, '.mat')
    oedb.datasetnames{u} = ff;
else
    oedb.datasetnames{u} = [ff ee];
end

if ~isempty(whos('-file', datafile, 'oerec'))
    
    oedb.issimulation(u) = false;
    if numel(oedb.datasetnames{u}) > 8 && strcmpi(oedb.datasetnames{u}(1:8),'dataset_')
        oedb.datasetnames{u} = oedb.datasetnames{u}(9:end);
    end
    oedb.data{u} = oerec_load(oedb.sourcedatafiles{u}); %also updates version if necessary
    oedb.nneurons(u) = numel(oedb.data{u});
        
elseif ~isempty(whos('-file', datafile, 'oesim'))
    
    oedb.issimulation(u) = true;
    if numel(oedb.datasetnames{u}) > 11 && strcmpi(oedb.datasetnames{u}(1:11),'simulation_')
        oedb.datasetnames{u} = oedb.datasetnames{u}(12:end);
    end
    oedb.data{u} = oesim_load(oedb.sourcedatafiles{u}); %also updates version if necessary
    oedb.nneurons(u) = numel(oedb.data{u}.oerec);
    
else
    
    error('invalid file');
    
end
assert(~isempty(oedb.data{u}), ['failed to load data from file ' oedb.sourcedatafiles{u}]);

oedb.expandneuron{u} = true(1, oedb.nneurons(u));
[oedb.neuronnames{u}, oedb.segmentnames{u}, oedb.segcode{u}, oedb.neuroncode{u}] = deal(cell(1, oedb.nneurons(u)));
oedb.nsegments{u} = nan(1, oedb.nneurons(u));
for v = 1:oedb.nneurons(u)
    
    nneuronstotal = nneuronstotal + 1;
    
    if oedb.issimulation(u)
        oerec = oedb.data{u}.oerec(v);
        oesim = oedb.data{u}; oesim.oerec = oesim.oerec(v);
    else
        oerec = oedb.data{u}(v);
    end
    oedb.nsegments{u}(v) = numel(oerec.data);
    if oedb.nsegments{u}(v) < 2
        oedb.expandneuron{u}(v) = false;
    end
    
    oedb.neuronnames{u}{v} = oedb_generate_neuron_name(oerec, v);
    oedb.neuroncode{u}{v} = num2str(nneuronstotal, 9);
    if ~oedb.opt.storedatainmemory
        if oedb.issimulation(u)
            save([fileinfo.tmpdir 's' oedb.neuroncode{u}{v} '.mat'],'oesim');
        else
            save([fileinfo.tmpdir 'n' oedb.neuroncode{u}{v} '.mat'],'oerec');
        end
    end
    
    ifiles = {oerec.data.imagefilepartialpath};
    ifiles_fullpath = {oerec.data.imagefilepartialpath};
    ifiles(cellfun(@isempty, ifiles)) = ifiles_fullpath(cellfun(@isempty, ifiles));
    for w = 1:oedb.nsegments{u}(v)
        if ~isempty(ifiles{w})
            [~,ifiles{w},~] = fileparts(ifiles{w});
        end
    end
    [~,~,ifileind] = unique(ifiles);
    
    [oedb.segcode{u}{v}, oedb.segmentnames{u}{v}] = deal(cell(1, oedb.nsegments{u}(v)));
    for w = 1:oedb.nsegments{u}(v)
        nsegtotal = nsegtotal + 1;
        oedb.segcode{u}{v}{w} = num2str(nsegtotal, 9);
        if ~isempty(oerec.data(w).imagefile)
            oedb.segmentnames{u}{v}{w} = ifiles{w};
            if sum(ifileind == ifileind(w)) > 1
                ifilesegind = find(w == find(ifileind == ifileind(w)));
                oedb.segmentnames{u}{v}{w} = [oedb.segmentnames{u}{v}{w} ' (seg ' num2str(ifilesegind) ')'];
            end
        else
            oedb.segmentnames{u}{v}{w} = num2str(w);
        end
        nmatchingchars = strnmatchingchars(oedb.segmentnames{u}{v}{w}, oedb.neuronnames{u}{v});
        if nmatchingchars > 5
            oedb.segmentnames{u}{v}{w} = oedb.segmentnames{u}{v}{w}(nmatchingchars + 1:end); %remove neuron name from segment name
        end
        if numel(oedb.segmentnames{u}{v}{w}) > 1 && ismember(oedb.segmentnames{u}{v}{w}(1), '/\_:., ')
            oedb.segmentnames{u}{v}{w} = oedb.segmentnames{u}{v}{w}(2:end);
        end
    end
    assert(nsegtotal < oedb_maxtotalsegs, ['oedatabrowser is restricted to less than ' num2str(oedb_maxtotalsegs) ' data segments']);
    
end
if ~oedb.opt.storedatainmemory
    
    oedb.data{u} = [];
    
end
oedb.ndatasets = u;
oedb.expanddataset(u) = true;
oedb.samedatasize = samedatasizemat(oedb);  % keep track of whether each pair of datasets can be compared fluorscence value by fluorescence value