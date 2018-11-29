function [r, alginfo] = apdet_constrainedfoopsi(oerec, params, opts)
%[results, alginfo] = apdet_constrainedfoopsi(oerec, params, opts)
%
%https://github.com/epnev/constrained-foopsi
assert(nargin == 3, '3 inputs required');

assert(exist('constrained_foopsi.m', 'file') == 2, 'constrained_foopsi.m was not found on the Matlab path. see https://github.com/epnev/constrained-foopsi');

nrecs = numel(oerec);
ndatasegments = arrayfun(@(o) numel(o.data), oerec);
r = apdet_resultsstruct(ndatasegments);

for n = 1:nrecs    
    
    for s = 1:ndatasegments(n)
        
        r{n}(s).params = params{n}{s};
        r{n}(s).opts = opts;
        
        [f, it, dt] = extract_data_fromoerec(oerec(n), s);
        f = reshape(f{1}, 1, []);
        it = it{1};
        
        [c,b,c1,g,sn,sp] = constrained_foopsi(f);
        
        r{n}(s).spikecount_times = it;
        r{n}(s).spikecounts = sp;
        r{n}(s).outputvars = struct('c',c,'b',b,'c1',c1,'g',g,'sn',sn);
        %g, sn, b and c1 could be made into parameters but for now let's not, since cfoopsi won't estimate them if they're passed. FIXME
        
    end
    
end

alginfo = empty_apdetalginfo;
alginfo.name = 'Constrained deconvolution';
alginfo.shortname = 'CFOOPSI';
alginfo.ismaximumlikelihood = true;
alginfo.isdeterministic = true;
alginfo.algorithm_version = 0.01;