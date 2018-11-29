function [fBL, blmask] = get_baseline_forfit(f, it, st, tau_s0, verbose)
%estimate baseline fluorescence
nseg = numel(f);
[fBL, blmask] = deal(cell(1, nseg));
for ii = 1:nseg
    %FIXME we should group together segments with the same session index when estimating baseline. maybe also force the same decay rate for the same neuron
    minboffset = 0;
    
    try
        
        [braw, ~, ~, blmask{ii}, boffset] = fit_unbinding_dynamics_withoffset(f{ii}, it{ii}, st{ii}, tau_s0, [], [], [], [], minboffset);
        fBL{ii} = braw + boffset;
        
    catch ex
        
        if verbose > 5
            
            fprintf('failed to estimate baseline for segment %d: %s\n', ii, ex.message);
            
        end
        
    end
    
end