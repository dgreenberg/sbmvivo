function [P, prevP, prevlik, conv] = estimate_parameters(V, P0, dt, F, nA2D)

points_per_segment = cellfun(@numel, F(:)');
skewness_each_segment = cellfun(@skewness, F(:)');
sd_each_segment = cellfun(@std, F(:)');
totalT = points_per_segment .* dt(:)';
init_F_win = 5; %sec
baselineguess_win = 120; %sec
baselineguess_percentile = 0.1;
init_F_npts = min(ceil(init_F_win ./ dt), points_per_segment);
baselineguess_npts = min(ceil(baselineguess_win ./ dt), points_per_segment);
[init_F, baselineguess] = deal(nan(1, numel(F)));
for u = 1:numel(F)
    
    init_F(u) = mean(F{u}(1:init_F_npts(u)));
    
    baselineguess_vals = sort(F{u}(1:baselineguess_npts(u)));
    percentile_index = ceil(baselineguess_npts(u) * baselineguess_percentile);
    baselineguess(u) = baselineguess_vals(percentile_index);
    
end
init_dff_guess = (init_F - baselineguess) ./ baselineguess;
shortsegs = totalT < max(totalT) / 2;
high_init_dff_segs = init_dff_guess > 0.2;

if sum(points_per_segment) > V.maxpoints_paramest
    
    %limit the number of time points used
    %use segments with greater skewness first
    [~, si] = sort(sd_each_segment, 'descend');
    
    %move short segments to the back of the list
    si = [si(~shortsegs(si)) si(shortsegs(si))];
    
    %move segments with high initial df/f values to the back of the list
    si = [si(~high_init_dff_segs(si)) si(high_init_dff_segs(si))];
    
    cs = cumsum(points_per_segment(si));
    segments_used = si([0 cs(1:end - 1)] < V.maxpoints_paramest);
    
    [Fpe, Ppe, dtpe, nA2Dpe] = deal(F(segments_used), P0(segments_used), dt(segments_used), nA2D(segments_used));
    
    if cs(numel(segments_used)) > V.maxpoints_paramest
        Fpe{end} = Fpe{end}(1:end - (cs(numel(segments_used)) - V.maxpoints_paramest));
    end
    
    if V.verbose > 2
        
        fprintf('Total data length of %d exceeds maximum length of %d, using partial data from %d segments to estimate neuron-specific parameters (all data will be used for AP inference).\n', ...
            sum(points_per_segment), V.maxpoints_paramest, numel(segments_used));
        
    end
    
    Vpe = V;
    fn = {'settingsindex', 'sessionindex', 'focalplaneindex'};
    for jj = 1:numel(fn)
        if isfield(V, fn{jj}) && ~isempty(V.(fn{jj}))
            Vpe.(fn{jj}) = V.(fn{jj})(segments_used);
        end
    end
    
else
    
    [Fpe, Ppe, dtpe, nA2Dpe, Vpe] = deal(F, P0, dt, nA2D, V);
    
end

switch V.parameterestimation
    
    case 'tm'
        
        [Ppe,prevPpe,prevlikpe] = est_params_temperedmode(Vpe, Ppe, dtpe, Fpe, nA2Dpe);
        
    otherwise
        
        error('unrecognized method');
        
end

if sum(points_per_segment) > V.maxpoints_paramest
    
    P = P0;
    prevP = repmat(P0, size(prevPpe, 1), 1);
    prevlik = repmat(empty_likstruct, size(prevPpe, 1), numel(P0));
    
    prevP(:,   segments_used) = prevPpe;
    prevlik(:, segments_used) = prevlikpe;
    
    for u = 1:numel(P)
        
        [P(u).R, P(u).S, P(u).fr] = deal(Ppe(1).R, Ppe(1).S, Ppe(1).fr);
        
    end
    
else
    
    [P, prevP, prevlik] = deal(Ppe, prevPpe, prevlikpe);
    
end

conv = false;