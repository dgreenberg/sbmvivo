function oerec = oerec_split_validF(oerec, minT)
expected_version = 0.01;
if ~exist('minT','var')
    minT = 10; %sec
end
oerec = update_oerec_version(oerec);
if numel(oerec) > 1
    for u = 1:numel(oerec)
        oerec(u) = oerec_split_validF(oerec(u), minT);
    end
    return;
elseif numel(oerec) == 0
    return;
end
assert(oerec.version == expected_version, 'oerec_split_validF needs to be updated');
newdata = oerec.data; newdata(:) = [];
for j = 1:numel(oerec.data)
    ok = oerec.data(j).f_mask & ~isnan(oerec.data(j).f);
    if ~any(ok)
        continue;
    end
    s = find(ok & [true; ~ok(1:end-1)]);
    e = find(ok & [~ok(2:end); true]);
    for k = 1:numel(s)
        
        data = oerecdata_extract_framerange(oerec.data(j), s(k):e(k));
        t = double(data.t) * data.t_scale_factor + data.t_offset;
        if max(t) - min(t) < minT
            
            [~, ff, ee] = fileparts(oerec.data(j).imagefile);
            warning('oerec_split_validF:segmenttooshort',['skipping segment ' num2str(k) ' for ' ff ee ', too short at ' num2str(max(data.t) - min(data.t)) ' sec']);
            continue;
            
        end
        
        data.notes{end+1} = ['segment ' num2str(k) ' / ' num2str(numel(s))];
        newdata = [newdata data]; %#ok<AGROW>
    end
end
oerec.data = newdata;