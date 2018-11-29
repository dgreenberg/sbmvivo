function [dt_eachsubseg, n_extrabins, bin_edges, t_states, stepsize_eachsubseg, nsteps_eachsubseg] = timinginfo( ...
    itsplit, target_substepsize, nsubseg, stsplit)
%determine step size, round each spike to the nearest modeling timestep, and get the time course for states/spikes
dt_eachsubseg = cellfun(@(v) diff(v(1:2)), itsplit); %for each subsegment

nsteps_eachsubseg = ceil(dt_eachsubseg / target_substepsize);
close_to_integer = abs(dt_eachsubseg / target_substepsize - round(dt_eachsubseg / target_substepsize)) < 0.01;
nsteps_eachsubseg(close_to_integer) = max(1, round(dt_eachsubseg(close_to_integer) / target_substepsize));

stepsize_eachsubseg = dt_eachsubseg ./ nsteps_eachsubseg;
[bin_edges, t_states] = deal(cell(1, nsubseg));
n_extrabins = zeros(1, nsubseg);
for u = 1:nsubseg
    
    bin_edges{u} = reshape( ...
        bsxfun(@plus, ...
        reshape(itsplit{u}, 1, []) - dt_eachsubseg(u), ...
        ((0:nsteps_eachsubseg(u) - 1)' - 0.5) * stepsize_eachsubseg(u)), ...
        [], 1);
    bin_edges{u}(end + 1, 1) = itsplit{u}(end) - stepsize_eachsubseg(u) / 2;
    
    %also include spikes before the first fluorescence measurement for this subsegment (e.g. the spike that begins the subsegment, or spikes before the beginning of a file)
    first_spiketime = stsplit{u}(1); %by definition there's at one least spike associated with each subsegment
    if first_spiketime < bin_edges{u}(1)
        
        n_extrabins(u) = ceil((bin_edges{u}(1) - first_spiketime) / stepsize_eachsubseg(u)) + 1;
        bin_edges{u} = [bin_edges{u}(1) + stepsize_eachsubseg(u) * (-n_extrabins(u):-1)'; bin_edges{u}];
        
    end
    
    t_states{u} = bin_edges{u}(1:end - 1) + 1.5 * stepsize_eachsubseg(u);
    
end