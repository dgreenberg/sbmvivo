function true_est_meanisi = oedb_burstinf_vs_n_and_meanisi(oedb, fileinfo, dsind, algind, maxisi, mintbefore, mintafter, maxtdiff, maxaps)
if ~exist('maxisi', 'var'), maxisi = 0.1; end
if ~exist('mintbefore', 'var'), mintbefore = 0.4; end %time before first AP in burst without any other APs
if ~exist('mintafter', 'var'), mintafter = 0.4; end %time after last AP in burst without any other APs
if ~exist('maxtdiff', 'var'), maxtdiff = 0.2; end % max time before first/after last AP in burst for inferred APs to be counted as detecting that AP
if ~exist('maxaps', 'var'), maxaps = 20; end % max APs in burst to be included in analysis
assert(maxisi < mintbefore && maxisi < mintafter && maxtdiff < mintbefore && maxtdiff < mintafter, 'maxisi and maxtdiff must be less than mintbefore and mintafter');

true_est_meanisi = cell(1, oedb.nneurons(dsind));
for neuronind = 1:oedb.nneurons(dsind)
    
    oerec = fetch_neurondata(oedb, fileinfo, dsind, neuronind);
    
    true_est_meanisi{neuronind} = zeros(0, 3);
    
    for segmentind = 1:oedb.nsegments{dsind}(neuronind)
        
        data = oerec.data(segmentind);
        results = fetch_results(oedb, fileinfo, dsind, neuronind, segmentind, algind);
        
        [spiketimes_true_for_stats, spiketimes_est_for_stats, weights_true_for_stats, weights_est_for_stats, twin_for_stats] = ...
            retrieve_spiketimes_for_stats(data, results);
        assert(all(weights_true_for_stats == 1), 'analysis of AP counts in bursts is not yet implemented for data with only binned AP counts available, but no AP times');
        
        st = spiketimes_true_for_stats(:);
        dst = diff(st);
        burst_starts = find([inf; dst(1:end - 1)] >= mintbefore & dst <= maxisi); %index of first AP in each burst we will include in this analysis
                
        nearedge = st(burst_starts) < twin_for_stats(1) + mintbefore | st(burst_starts) > twin_for_stats(2) - mintafter;
        burst_starts(nearedge) = [];
                
        next_true_est_meanisi = nan(numel(burst_starts), 3);
        valid = false(numel(burst_starts), 1);
        for j = 1:numel(burst_starts)
            
            burst_end = burst_starts(j) + find([dst(burst_starts(j):end); inf] > maxisi, 1) - 1; %index of last in this burst
            naps_true = burst_end - burst_starts(j) + 1;
            
            if st(burst_end) > twin_for_stats(2) - mintafter || ...
                    naps_true > maxaps || ...
                    (j < numel(burst_starts) && dst(burst_end) < mintafter) % too close to end of data or not enough time without APs after last AP in burst
                
                continue;
            
            end
            valid(j) = true;
            
            ii = spiketimes_est_for_stats >= st(burst_starts(j)) - maxtdiff & spiketimes_est_for_stats <= st(burst_end) + maxtdiff;
            naps_est = sum(weights_est_for_stats(ii));
            
            meanisi = (st(burst_end) - st(burst_starts(j))) / (naps_true - 1);
            
            assert(j == numel(burst_starts) || burst_starts(j + 1) > burst_end, 'invalid program state');
            
            next_true_est_meanisi(j, :) = [naps_true naps_est meanisi];                        
            
        end
        true_est_meanisi{neuronind} = [true_est_meanisi{neuronind}; next_true_est_meanisi(valid, :)];
        
    end
        
end