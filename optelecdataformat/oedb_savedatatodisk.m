function oedb_savedatatodisk(oedb, fileinfo)
nalgs = numel(oedb.algnames);
for u = 1:oedb.ndatasets
    for v = 1:oedb.nneurons(u)
        if oedb.issimulation(u)
            
            oesim = oedb.data{u}(v); %#ok<NASGU>
            save([fileinfo.tmpdir 's' oedb.neuroncode{u}{v} '.mat'],'oesim');
            
        else
            
            oerec = oedb.data{u}(v); %#ok<NASGU>
            save([fileinfo.tmpdir 'n' oedb.neuroncode{u}{v} '.mat'],'oerec');
            
        end
        for w = 1:oedb.nsegments{u}(v)
            for algind = 1:nalgs
                
                results = fetch_results(oedb, fileinfo, u, v, w, algind); %#ok<NASGU>
                rfile = rfilename(oedb, fileinfo, algind, u, v, w);
                save(rfile, 'results');
                
            end
        end
    end
end
