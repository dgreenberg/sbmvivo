function assign_defaults(varargin)
if mod(nargin,2)
    error('even number of inputs required for assign_defualts.m');
end

for k = 1:2:nargin
    if strcmp(varargin{k}(1),'\') %variable name preceded by a slash, indicating that we should evaluate the subsequent string in the base workspace, instead of assinging it
        varargin{k} = varargin{k}(2:end);
        if ~evalin('caller',['exist(''' varargin{k} ''', ''var'')'])
            if ~strcmp(varargin{k+1}(end),';')
                varargin{k+1}(end + 1) = ';';
            end
            evalin('caller', [varargin{k} ' = ' varargin{k+1}]);
        end
    else
        if ~evalin('caller',['exist(''' varargin{k} ''', ''var'')'])
            assignin('caller', varargin{k}, varargin{k+1});
        end
    end
end