function moments = empty_momentstruct(T, nsteps, V)
if ~exist('T','var')
    
    T = 0;
    
end
if ~exist('nsteps','var')
    
    nsteps = 1;
    
end

allstates = {'c', 's', 'x', 'n', 's0', 's1', 's2', 's3', 's4', 'b0', 'b1', 'r', 'expr', 'gp', 'Fpred', 'n'};

statenames = {'c', 'r', 'expr', 'gp', 'Fpred', 'n'};

if exist('V', 'var') && T > 0
    
    if strcmpi(V.model, '5s2b')
        
        statenames = [statenames {'s0', 's1', 's2', 's3', 's4', 'b0', 'b1'}];
        
    else
        
        error('unrecognized model: %s', V.model);
        
    end
    
end

moments = struct();
for k = 1:numel(allstates)
    
    sn = allstates{k};
    
    if ~ismember(sn, statenames)
        
        [moments.([sn '_mean']), moments.([sn '_sd'])] = deal(nan(1, 0));
        continue;
        
    end
    
    if ismember(sn, {'c', 's', 'x', 'n', 's0', 's1', 's2', 's3', 's4', 'b0', 'b1'})
        
        npts = T * nsteps;
        
    else
        
        npts = T;
        
    end
    [moments.([sn '_mean']), moments.([sn '_sd'])] = deal(nan(1, npts));
    
end
%moments.sexprsq = deal(nan(1, T)); %sufficient statistic need for estimation of R
moments.nsteps = nsteps;
moments = orderfields(moments);