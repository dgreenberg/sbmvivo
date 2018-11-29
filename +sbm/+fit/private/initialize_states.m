function [states, states_eq] = initialize_states(P, t_states, modelname)
nsubseg = numel(t_states);
if isfield(P, 'Btot')
    
    nbuffers = arrayfun(@(v) numel(v.Btot), P);
    assert(all(nbuffers == nbuffers(1)), 'all segments must have the same number of buffers');
    nbuffers = nbuffers(1);
    
else
    
    nbuffers = 0;
    
end
if modelname(1) == '5'
    
    nstates = 6 + nbuffers;  % calcium, macromolecule binding states, calcium-bound buffer states
    
else
    
    nstates = 3 + nbuffers; % calcium, intermediate and saturated macromolecule binding states, calcium-bound buffer states
    
end

[states] = deal(cell(1, nsubseg));
states_eq = nan(nstates, nsubseg);

for k = 1:nsubseg
    
    % note we don't use deal here, causes problems for mex access later!!
    states{k} = nan(nstates, numel(t_states{k}));
    
end