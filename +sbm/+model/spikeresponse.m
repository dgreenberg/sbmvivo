function [brightness, states, tt, spike_time] = spikeresponse(P,opts,dt,T,ns)
assert(numel(P) == 1, 'one set of parameters only');
if ~exist('ns','var') || isempty(ns)
    ns = 1;
end
if ~exist('dt','var') || isempty(dt)
    dt = 0.005;
end
if ~exist('T','var') || isempty(T)
    T = 2 + ceil(5 / dt); %fixme bad notation
end

if isfield(opts, 'substepsize')
    
    nsteps = ceil(dt / opts.substepsize); %10 ms increments or so
    
else
    
    nsteps = 1;
    
end

nsvec = zeros(T * nsteps, 1);
nsvec(nsteps + 1) = ns;

[brightness, states, ~, tt] = sbm.model.spiketrainresponse(P,opts,dt,nsvec);

spike_time = tt(nsteps);  % spike times rae shifted back one step relative to state times