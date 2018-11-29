function M = gfitspiketrain(moments, V, dt)
if numel(moments) > 1
    for k = 1:numel(moments)
        M(k) = gfitspiketrain(moments(k), V, dt(k)); %#ok<AGROW>
    end
    return;
end

stepsize = dt / moments.nsteps;
totalT = numel(moments.n_mean) * stepsize;
estfr = sum(moments.n_mean) / totalT;

tt = stepsize * (1:numel(moments.n_mean)); %times for states (not spikes!) relative to one dt before first fluorescence measurement
if estfr > V.gfitst_maxfr && sum(moments.n_mean) > 40 %too many spikes, so fitting would take too long
    
    warning('not estimating spike train for data with inferred firing rate of %f', estfr);
    [st, stsd] = deal([]);
    gfit = nan(size(moments.n_mean));
    
else
    
    maxsd = (V.maxjitter_ms / 1000) / stepsize;
    [st, stsd, gfit] = moments2jitter(moments.n_mean, maxsd);
    
end
M.gst = orderfields(struct('mu',st * stepsize, 'sd', stsd * stepsize, 'gfit', gfit, 'tt', tt));
M.version = sbm.version;