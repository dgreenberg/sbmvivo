function [aptimes, aptimes_window, params, opts] = infer(f, it, indicator, params, opts, nA2D)
% simplest call:
% aptimes = infer(f, it, indicator)
% full call:
% [aptimes, aptimes_window, params, opts] = infer(f, dt, params, opts, nA2D)
%
% REQUIRED INPUTS:
% f -- cell array of fluorescence values, or a vector for a single segment
% it -- image times, same structure as f
% indicator -- string, implemented for gcamp6s so far
%
% OPTIONAL INPUTS
% params -- strucutre of parameter values inferred from a previous run or known by other means. Can be left empty.
%           params can be a structure array the same size as F when dealing with multiple data segments
% opts -- structure array of options.
%         for GPU computation, set opts.usegpu = true
%
% OUTPUTS
% aptimes -- cell array of inferred AP times for each data segment
% aptimes_window -- time window for which inference was carried out
% params, opts, -- structures containing the parameters and opts used for inference.
%                  these will include both user-supplied values, default values,
%                  and values that were assigned automatically based on the data
%
% This function is a wrapper for sbm.run_alg, which provides more fine
% grained acces to inputs/outputs.

if isnumeric(f)
    f = {f};
    assert(isnumeric(it));
    it = {it};
end
nsegments = numel(f);
assert(numel(it) == nsegments && all(cellfun(@numel, it) == cellfun(@numel, f)), 'f and it do not match');

if ~exist('indicator', 'var'), indicator = ''; end
if ~exist('params', 'var'),    params = []; end
if ~exist('opts', 'var'),      opts = []; end
if ~exist('nA2D', 'var'),      nA2D = []; end

dt = cellfun(@(x) median(diff(x)), it);
max_deviation = cellfun(@(x) max(abs(diff(x) - median(diff(x)))) / median(diff(x)), it);
if any(max_deviation > 1.02)
    
    warning('Some fluorescence measurements are not evenly spaced in time; max asbsolute deviation from median %g %%', 100 * max(max_deviation));
    
end

[M, params, opts, ~, ~, ~, moments] = sbm.run_alg(f, dt, indicator, opts, params, [], [], nA2D);
aptimes = cell(1, nsegments);
aptimes_window = nan(nsegments, 2);
for s = 1:nsegments
    
    nframes = numel(it{s});
    nspikes = numel(M(s).gst.mu);
    st = nan(nspikes, 1);
    
    %the states and observations have to be shifted back by stepsize to get spike times:
    stepsize = dt(s) / moments(s).nsteps;
    st_steps = M(s).gst.mu / stepsize - 1;  % get the times of non-discretized detected APs in simulation steps. zero is exactly one frame (dt) before the first fluorescence measurement
    st_frames = st_steps / moments(s).nsteps;  % get the times of non-discretized detected APs in frames. zero is exactly one frame (dt) before the first fluorescence measurement
    st_frames_int = floor(st_frames);
    st_frames_frac = st_frames - st_frames_int;
    
    before_first_measurement = st_frames < 1;
    st(before_first_measurement) = it{s}(1) + (st_frames(before_first_measurement) - 1) * dt(s);
    during_measurement = st_frames >= 1 & st_frames <= nframes;
    st(during_measurement) = it{s}(st_frames_int(during_measurement)) + st_frames_frac(during_measurement) * dt(s);
    
    after_last_measurement = st_frames > nframes;
    st(after_last_measurement) = it{s}(end) + (st_frames(after_last_measurement) - nframes) * dt(s);
    
    aptimes{s} = st;
    aptimes_window(s, :) = [it{s}(1) - dt(s), it{s}(end) - stepsize] + 0.5 * stepsize * [-1 1];
    
end