function oedb = oedb_loaddatafromdisk(oedb, fileinfo)
nalgs = numel(oedb.algnames);
for u = 1:oedb.ndatasets
    for n = 1:oedb.nneurons(u)
        if oedb.issimulation(u)
            
            oedb.data{u}(n) = oesim_load([fileinfo.tmpdir 's' oedb.neuroncode{u}{n} '.mat']); %also updates version if necessary
            
        else
            
            oedb.data{u}(n) = oerec_load([fileinfo.tmpdir 'n' oedb.neuroncode{u}{n} '.mat']); %also updates version if necessary
            
        end
        for s = 1:oedb.nsegments{u}(n)
            for algind = 1:nalgs
                
                if exist([fileinfo.tmpdir 'r' num2str(algind) '_' oedb.segcode{u}{n}{s} '.mat'], 'file')
                    
                    load([fileinfo.tmpdir 'r' num2str(algind) '_' oedb.segcode{u}{n}{s} '.mat'],'results');
                    
                else
                    
                    results = apdet_resultsstruct(1); %initialize results struct for one datasegment
                    results = results{1};
                    
                end
                
                oedb.results{algind, u}{n}(s) = results;
                
            end
        end
    end
end