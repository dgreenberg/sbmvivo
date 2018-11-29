function canonlinfit_save_figures(Pout, segmentind_split, neuronnames_eachsegment, states_combo, states_eq_eachsubseg, t_states_combo, ...
    it, fBL, f, fhat_combo, blmask_combo, rms_combo, R, Rind, L, pi, statenames, perispikewin, basefigsavedir)
nseg = numel(f);
if isfield(Pout(1), 'concentration_scale_to_nM')
    
    states_scalefac = Pout(1).concentration_scale_to_nM;
    
else
    
    states_scalefac = 1;
    
end

fignames = cell(1, nseg);
for u = 1:nseg
    
    if ~any(segmentind_split == u), continue; end %segment had no spikes or no baseline, and wasn't used
    
    fignames{u} = [neuronnames_eachsegment{u} ' segment ' num2str(u)];
    
    plot_canonlinfit(states_combo{u}, states_eq_eachsubseg(:, find(segmentind_split == u, 1)), t_states_combo{u}, ...
        it{u}, fBL{u}, f{u}, fhat_combo{u}, fignames{u}, blmask_combo{u}, ...
        rms_combo(u), sum(segmentind_split == u), 0, R(Rind(u)), L(pi.S(u)), ...
        statenames, perispikewin, basefigsavedir, states_scalefac);
    
end

if nseg > 1
    esegfig = figure('windowbuttondownfcn', @(figh, eventdata) esegcombofib_buttondowncb(figh, eventdata, fignames));
    plot(rms_combo, 'parent', axes('parent', esegfig));
end