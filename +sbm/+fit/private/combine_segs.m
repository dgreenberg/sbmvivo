function [fBL_combo, blmask_combo, fhat_combo, states_combo, t_states_combo, eseg_combo] = combine_segs( ...
    f, fBL_split, blmask_split, fhat_split, states_split, segmentind_split, timeind, t_states, eseg_split)
nseg = numel(f);
nseg_split = numel(fBL_split);
nstates = size(states_split{1}, 1);

[fBL_combo, blmask_combo, fhat_combo, states_combo, t_states_combo] = deal(cell(1, nseg));
eseg_combo = zeros(1, nseg);
for v = 1:nseg
    
    [fBL_combo{v}, fhat_combo{v},  states_combo{v}] = deal(nan(size(f{v})));
    blmask_combo{v} = false(size(f{v}));
    states_combo{v} = zeros(nstates, 0);
    t_states_combo{v} = [];
    
    if ~any(segmentind_split == v)
        
        eseg_combo(v) = nan;
        
    end
    
end
for u = 1:nseg_split
    
    v = segmentind_split(u);
    
    fBL_combo{v}(timeind{u})    = fBL_split{u};
    fhat_combo{v}(timeind{u})   = fhat_split{u};
    blmask_combo{v}(timeind{u}) = blmask_split{u};
    states_combo{v} = [states_combo{v} states_split{u} nan(nstates, 1)];
    t_states_combo{v} = [t_states_combo{v}; t_states{u}; nan];
    eseg_combo(v) = eseg_combo(v) + eseg_split(u);
    
end