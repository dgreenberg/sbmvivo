function y = copy_matlab_var(x)
% y = copy_matlab_var(x)
%
% Creates a copy with memory in a different location, so mex files
% operating on the copy will not affect the data in the original.
% 
% Supports numeric, char and logical arrays, as well as cells and structs
% with any type of nesting. Other classes are not supported.
if isnumeric(x) || isa(x, 'char') || isa(x, 'logical')
    
    y = x;
    y(1) = y(1);
    
elseif isa(x, 'struct')
    
    fn = fieldnames(x);
    y = struct();
    for k = 1:numel(fn)
       
        y.(fn{k}) = copy_matlab_var(x.(fn{k}));
        
    end
    
elseif isa(x, 'cell')

    y = cell(size(x));
    for k = 1:numel(x)
        
        y{k} = copy_matlab_var(x{k});
        
    end
    
else
        
    error('copy_matlab_var not yet implemented for class: %s', class(x));
    
end