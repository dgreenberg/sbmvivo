function [f_split,it_split,st_split,fBL_split,timeind,segmentind_split,blmask_split] = split_segs(f, it, st, fBL, blmask, perispikewin)
minsegt = 0.3; %sec
nseg = numel(f);
[f_split, it_split, st_split, fBL_split, timeind, blmask_split] = deal({});
[segmentind_split] = deal([]);
for u = 1:nseg
    if isempty(blmask{u}), continue; end %failed to get baseline
    ok = false(size(it{u}));
    for j = 1:numel(st{u})
        ok(it{u} >= st{u}(j) + perispikewin(1) & it{u} <= st{u}(j) + perispikewin(2)) = true;
    end
    s = find(ok & [true; ~ok(1:end-1)]);
    e = find(ok & [~ok(2:end); true]);
    tooshort = it{u}(e) - it{u}(s) < minsegt;
    s(tooshort) = []; e(tooshort) = [];
    for k = 1:numel(s)
        f_split{1, end + 1} = reshape(f{u}(s(k):e(k)), [], 1); %#ok<*AGROW>
        it_split{1, end + 1} = it{u}(s(k):e(k));
        st_split{1, end + 1} = st{u}(st{u} >= it{u}(s(k)) - perispikewin(2) & st{u} <= it{u}(e(k)) + perispikewin(1));        
        segmentind_split(1, end + 1) = u;
        timeind{1, end + 1} = s(k):e(k);
        fBL_split{1, end + 1} = fBL{u}(s(k):e(k));
        blmask_split{1, end + 1} = blmask{u}(s(k):e(k));
    end
end