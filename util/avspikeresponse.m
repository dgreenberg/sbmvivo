function [av, tav, av_sd, av_sem, responsedata, naps_isolated, naps_total, fr_total, T_total, tevent, segind, stdiff] = ...
    avspikeresponse(oerec, mintbefore, extractdt, subtractprev, mintbefore_singlespike, subwin)
if ~exist('subtractprev','var') || isempty(subtractprev)
    subtractprev = true;
end
if ~exist('mintbefore','var') || isempty(mintbefore)
    mintbefore = 1.5;
end
if ~exist('mintbefore_singlespike','var') || isempty(mintbefore_singlespike)
    mintbefore_singlespike = min(mintbefore, max(1, mintbefore / 2));
end
if ~exist('extractdt','var') || isempty(extractdt)
    extractdt = 0.075;
end
if ~exist('subwin', 'var') || isempty(subwin)
    subwin = [-.2 -.001];
end
if numel(oerec) > 1
    
    if numel(mintbefore) == 1, mintbefore = repmat(1, mintbefore, numel(oerec)); end
    if numel(mintbefore_singlespike) == 1, mintbefore_singlespike = repmat(1, mintbefore_singlespike, numel(oerec)); end
    
    [av, tav, av_sd, av_sem, responsedata, naps_isolated, tevent, segind, stdiff] = deal(cell(1, numel(oerec)));
    [naps_total, fr_total, T_total] = deal(nan(1, numel(oerec)));
    for j = 1:numel(oerec)
        
        [av{j}, tav{j}, av_sd{j}, av_sem{j}, responsedata{j}, naps_isolated{j}, naps_total(j), fr_total(j), T_total(j), tevent{j}, segind{j}, stdiff{j}] = ...
            avspikeresponse(oerec(j), mintbefore(j), extractdt, subtractprev, mintbefore_singlespike(j));
        
    end
    return;
    
end

minboffset = 0;
tau_s0 = 0.5; %sec
mintafter = 0.15; %sec
isolation_win = 0.2;
extractwin = [-1 4];
subwin(1) = max(subwin(1), -mintbefore_singlespike);
avedges = extractdt * (floor(extractwin(1) / extractdt):ceil(extractwin(2) / extractdt))';
tav = avedges(1:end-1)' + extractdt / 2;
naps_total = 0; T_total = 0;
[responsedata,tevent,segind,stdiff] = deal({});
nspikes = [];

for seg = 1:numel(oerec.data)
    
    f = reshape(double(oerec.data(seg).f) * oerec.data(seg).f_scale_factor + oerec.data(seg).f_offset,[],1);
    t = reshape(double(oerec.data(seg).t) * oerec.data(seg).t_scale_factor + oerec.data(seg).t_offset,[],1);
    dt = mean(diff(t));
    T_total = T_total + max(t) - min(t) + dt;
    if oerec.data(seg).spiketimes_present
        
        st = reshape(oerec.data(seg).spiketimes,[],1);
        st = sort(st);
        st_used = st(st >= t(1) - dt & st <= t(end)); % included in the average. we still need to keep track of spikes before the first fluorescence value        
        
    elseif oerec.data(seg).apcounts_present
        
        st = [];
        for f = reshape(find(oerec.data(seg).apcounts), 1, [])
            ns = oerec.data(seg).apcounts(f);
            if f == 1
                tprev = t(1) - dt;
            else
                tprev = t(f - 1);
            end
            st = [st; tprev + dt * (1:ns)' / (ns + 1)]; %#ok<AGROW> %space out spikes evenly across frame when we don't know their exact times
        end
        st_used = st;
        
    else
        
        warning('ephys must be available, skipping segment %d / %d', seg, numel(oerec.data));
        continue;
        
    end
    assert(numel(st) == numel(unique(st)), 'duplicate spike time');
    
    naps_total = naps_total + numel(st_used);  % don't include spikes from excluded data, since this may include e.g. the cell dying etc.
    try
        
        [braw, tau_s, s, ~, boffset] = fit_unbinding_dynamics_withoffset(f, t, st, tau_s0, [], [], [], [], minboffset, 0);  % keep all spike times in mind
        b = braw + boffset;        
        
    catch ex
        
        [~, ff, ee] = fileparts(oerec.data(seg).imagefile);
        segname = [ff ee];
        if ~isempty(oerec.data(seg).notes) && numel(oerec.data) > 1
            segname = [segname '(' oerec.data(seg).notes{1} ')']; %#ok<AGROW>
        end
        warning('avspikeresponse:baselinefailed',['Failed to calculate baseline for ' segname ':' ex.message]);
        continue;
        
    end
    dff = (f - b) ./ b;
    
    for j = 1:numel(st_used)
        
        nextst = st_used(j);
        
        ost = st(st ~= nextst);  % other spike times, including those outside of imaging period
        preceding_spikes = ost(ost <= nextst & ost >= nextst - mintbefore);
        if ~isempty(preceding_spikes)
            if numel(preceding_spikes) > 1 || any(nextst - preceding_spikes <= mintbefore_singlespike)
                continue;
            end
        end
        if any(ost > nextst + isolation_win & ost < nextst + mintafter)
            continue;
        end
        
        if subtractprev && ~isnan(tau_s)
            
            valid_range_prevsind = find(t < nextst & t > max([-inf; st(st < nextst)]));
            prevsind = valid_range_prevsind(find(~isnan(s(valid_range_prevsind)), 1, 'last'));
            if isempty(prevsind)
                continue;
            end
            
        end
        
        spikes_in_burst = [nextst; ost(ost >= nextst & ost <= nextst + isolation_win)];
        
        nspikes(1, end + 1) = numel(spikes_in_burst); %#ok<AGROW>
        nextextractwin = extractwin + nextst;
        if any(ost > nextst + isolation_win) %only include time points before the next spike not part of the current burst
            
            nextextractwin(2) = min(nextextractwin(2), min(ost(ost > nextst + isolation_win)));
            
        end
        ii = t >= nextextractwin(1) & t <= nextextractwin(2);
        if numel(responsedata) < nspikes(1, end)
            
            [responsedata{nspikes(1, end)}, tevent{nspikes(1, end)}, segind{nspikes(1, end)}, stdiff{nspikes(1, end)}] = deal({});
            
        end
        if subtractprev && ~isnan(tau_s)
            
            fii = find(ii);
            expected_s_ii = s(ii);
            expected_s_ii(fii > prevsind) = s(prevsind) * exp(-dt * (fii(fii>prevsind) - prevsind) / tau_s); %expected s for frames in range described by ii given previous activity only
            %note that far enough back, expected_s_ii and thus nextdff may be nan due to nonzero calcium after a spike, which we do not try to model here
            nextdff = dff(ii) - expected_s_ii;
            
        else
            
            nextdff = dff(ii);
            
        end
        subval = mean_nonnan(nextdff(t(ii) >= nextst + subwin(1) & t(ii) <= nextst + subwin(2)));
        if ~isnan(subval)
            
            nextdff = nextdff - subval;
            
        end
        
        responsedata{nspikes(1, end)}{1, end + 1} = [t(ii) - nextst nextdff];
        tevent{nspikes(1, end)}{1, end + 1} = nextst;
        segind{nspikes(1,end)}{1, end + 1} = seg;
        stdiff{nspikes(1,end)}{1, end + 1} = diff(spikes_in_burst);
        
    end
    
end
naps_isolated = histc(nspikes, 1:max(nspikes));
fr_total = naps_total / T_total;
[av, av_sd, av_sem] = deal(cell(1, max(nspikes)));
for j = 1:max(nspikes)
    if isempty(responsedata{j})
        responsedata{j} = {};
        [av{j},av_sd{j},av_sem{j}] = deal(nan(numel(tav), 1));
        continue;
    end
    responsedatamat = cat(1, responsedata{j}{:});
    [av{j},av_sd{j},av_sem{j}] = histfunc(responsedatamat(:,1), responsedatamat(:,2),avedges);
    av{j} = av{j}(1:end-1); av_sd{j} = av_sd{j}(1:end-1); av_sem{j} = av_sem{j}(1:end-1);
end