function [r, alginfo] = apdet_foopsi(oerec, params, opts)
%[results, alginfo] = apdet_foopsi(oerec, params, opts)
%
%Fast nonnegative
%oerec is an array of optical / electrical datasets
%Vogelstein et al. 2010
assert(nargin == 3, '3 inputs required');

nrecs = numel(oerec);
ndatasegments = arrayfun(@(o) numel(o.data), oerec);
r = apdet_resultsstruct(ndatasegments);

for n = 1:nrecs
    
    for s = 1:ndatasegments(n)
        
        r{n}(s).params = params{n}{s};
        r{n}(s).opts = opts;
        
        [f,it,dt] = extract_data_fromoerec(oerec, s);
        f = reshape(f{1}, 1, []);
        it = it{1};
        
        assert(~any(isnan(f)),'not implemented for data with NaNs');
        T = numel(f);
        
        %detrend and normalize data
        f_fast = z1(detrend(f));
        
        opts_forfast = orderfields(struct( ...
            'fast_plot',    0 ...
            ,'dt',           dt ...
            ,'T',            T ...
            ,'est_sig',      1 ...
            ,'est_b',        0 ...
            ));
        
        %assign default parameters if missing
        if ~isfield(r{n}(s).params, 'gam') || ~isfield(r{n}(s).params, 'a')
            
            tau_c   = 0.5;
            r{n}(s).params.gam = 1 - dt / tau_c;
            r{n}(s).params.a   = dt / tau_c;
            
        end
        
        r{n}(s).spikecounts = fast_oopsi(f_fast, opts_forfast, r{n}(s).params);
        r{n}(s).spikecount_times = it;
        
    end
end

alginfo = empty_apdetalginfo;
alginfo.name = 'Fast Nonnegative Deconvolution';
alginfo.shortname = 'FOOPSI';
alginfo.ismaximumlikelihood = true;
alginfo.isdeterministic = true;
alginfo.algorithm_version = 0.01;


function x = z1(y)
% linear normalize between 0 and 1
x = (y-min(y(:)))/(max(y(:))-min(y(:)))+eps;