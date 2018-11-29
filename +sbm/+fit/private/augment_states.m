function [states_a, states_eq_a] = augment_states(states, states_eq, P_eachsubseg, modelname)
nseg = numel(states);
nbuffers = arrayfun(@(v) numel(v.Btot), P_eachsubseg);
assert(all(nbuffers == nbuffers(1)), 'all segments must have the same number of buffers');
nbuffers = nbuffers(1);

states_a = cell(size(states));

if modelname(1) == '5'
    
    states_eq_a = nan(size(states_eq, 1) + nbuffers, nseg);
    
else
    
    states_eq_a = nan(size(states_eq, 1) + 1 + nbuffers, nseg);
    
end

for jj = 1:nseg
    
    T = size(states{jj}, 2);
    S = P_eachsubseg(jj).S;
    
    if nbuffers > 0
        
        Btot = P_eachsubseg(jj).Btot(:);
        CB = states{jj}(end - nbuffers + 1:end, :); %concentration of calcium-bound state of each buffer at each time step
        Bfree = bsxfun(@minus, Btot, CB); %concentration of calcium-free state of each buffer at each time step
        bufferstates = reshape([Bfree CB]', T, [])'; %1st row calcium free buffer 1, 2nd row calcium bound buffer 1, 3rd row calcium free buffer 2, etc.
        
        CB_eq = states_eq(end - nbuffers + 1:end, jj);
        bufferstates_eq = reshape([Btot - CB_eq CB_eq]', [], 1);
        
    else
        
        bufferstates = zeros(0, T);
        bufferstates_eq = zeros(0, 1);
        
    end
    
    if modelname(1) == '5'
        
        mm_states = states{jj}(2:end - nbuffers, :);
        mm_states_eq = states_eq(2:end - nbuffers, jj);
        
    else
        
        mm_states = [S - sum(states{jj}(2:3, :), 1); states{jj}(2:3, :)];
        mm_states_eq = [S - sum(states_eq(2:3, jj), 1); states_eq(2:3, jj)];
        
    end
    
    states_a{jj} = [states{jj}(1, :); ... %free calcium
        mm_states;
        bufferstates];
    
    states_eq_a(:, jj) = [states_eq(1, jj); ... %free calcium
        mm_states_eq; ... %calcium-free indicator
        bufferstates_eq];
    
end