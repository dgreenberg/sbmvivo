function v = checkvectorinput(v, nelements, name)
if ~exist('name','var') || isempty('name')
    name = 'input';
end
assert(isnumeric(v) && numel(v) == nelements && ~any(isinf(v) | (imag(v) ~= 0)), '%s must be a vector of non-complex, non-inf doubles', name);
v = double(reshape(v, 1, nelements));