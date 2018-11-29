function opts = assignopts(opts)
if ~isfield(opts, 'modelname')
    
    opts.modelname = '5sb';
    
end
dopts = defaultopts(opts.modelname);
fn = fieldnames(dopts);
for k = 1:numel(fn)
    
    if ~isfield(opts, fn{k})
        
        opts.(fn{k}) = dopts.(fn{k});
        
    end
    
end
%fast and slow buffers
if ~isfield(opts, 'min_tau_B')
    
    opts.min_tau_B = 1e-3 * ones(1, opts.nbuffers);
    
end
if ~isfield(opts, 'max_tau_B')
    
    opts.max_tau_B = 30 * ones(1, opts.nbuffers);
    
end