function varargout = tps_mlspikes_wrapper(varargin) %#ok<STOUT>
mlspikedir = fileparts(which('tps_mlspikes.m'));
brickdir = [mlspikedir filesep '..' filesep 'brick'];
addpath(brickdir);
try
    
    if nargout > 0
        
        argsoutstr = '[';
    
    else
        
        argsoutstr = '';
        
    end
        
    for k = 1:nargout
        
        argsoutstr = [argsoutstr 'varargout{' num2str(k) '} ']; %#ok<AGROW>
        
    end
    if nargout > 0
        
        argsoutstr = [argsoutstr '] = '];
        
    end
    
    argsinstr = '';
    for k = 1:nargin
        
        if k > 1
            
            argsinstr = [argsinstr ',']; %#ok<AGROW>
            
        end
        
        argsinstr = [argsinstr ' varargin{' num2str(k) '}']; %#ok<AGROW>
        
    end
    
    eval([argsoutstr 'tps_mlspikes(' argsinstr ');']);
    rmpath(brickdir);
    
catch ex
    
    rmpath(brickdir);
    rethrow(ex);
    
end