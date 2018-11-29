function [r, alginfo] = apdet_dpncls(oerec, params, opts)
%[results, alginfo] = apdet_dpncls(oerec, params, opts)
%
%dynamic programming / nonconvex least squares
%oerec is an array of optical / electrical datasets
% %Greenberg et al. 2008
assert(nargin == 3, '3 inputs required');

nrecs = numel(oerec);
ndatasegments = arrayfun(@(o) numel(o.data), oerec);
r = apdet_resultsstruct(ndatasegments);

for n = 1:nrecs
    
    for s = 1:ndatasegments(n)
                
        r{n}(s).params = params{n}{s};
        r{n}(s).opts = opts;
        
        [f,it,dt,~,indicatorstring] = extract_data_fromoerec(oerec(n), s);
        f = reshape(f{1}, [], 1);
        it = it{1};
        indicatorstring = indicatorstring{1};
        
        %assign parameters if not supplied as an input
        if isempty(r{n}(s).params) || isempty(fieldnames(r{n}(s).params)) %fixme merge structs
            
            switch indicatorstring
                case 'ogb1-am'
                    
                    r{n}(s).params = dpncls_default_spikefind_params('OGB1-AM'); %FIXME make these params
                    
                otherwise
                    
                    warning('dpncls:indicatornotrecognized','unrecognized indicator, using parameters for ogb1-am');
                    r{n}(s).params = dpncls_default_spikefind_params('OGB1-AM');
                    
            end
            
        end
        r{n}(s).outputvars.blwin = 10; r{n}(s).outputvars.nearnoregthresh = 0.17; r{n}(s).outputvars.afterspikewin = 1.75;
        
        [amp, spmodel, finalbl] = dpncls.spikemodel_paramwrapper(f, dt, r{n}(s).params, inf);
        
        r{n}(s).spikecounts = round(amp / 9.75);
        r{n}(s).spikecount_times = it;
        r{n}(s).outputvars.finalbl = finalbl;
        
    end
end

alginfo = empty_apdetalginfo;
alginfo.name = 'Dynamic programming / nonconvex least squares';
alginfo.shortname = 'DPNCLS';
alginfo.ismaximumlikelihood = true;
alginfo.isdeterministic = true;
alginfo.algorithm_version = 0.011;

function params = dpncls_default_spikefind_params(type)
assign_defaults('type','OGB1-AM')
switch type
    case 'OGB1-AM'
        params = struct(...
            'regwin',                   [-1 3],...
            'fitwin',                   1.75,...
            'minth',                    10,...
            'minindamp',                8,...
            'minpeakderiv',             2,...
            'minpeakdff',               6.1,...
            'tau',                      0.5,...
            'blwin',                    20,...
            'blgsd',                    0.3,...
            'minregpeak',               7.45,...
            'pval',                     50,...
            'nearnoregthresh',          0.27,...
            'highthresh',               9,...
            'minpeakderivhigh',         2,...
            'minstradle',               -3.15,...
            'mindepbestamp',            6.6,...
            'reblsd',                   0.75,...
            'reblsd2',                  2,...
            'reblwin',                  20,...
            'mindepmean',               6.55,...
            'minpeakdoubderiv',         0.8,...
            'mindepdoubderiv',          3,...
            'areawin',                  0.4,...
            'dffperspike',              9.75,...
            'version',                  1.0,...
            'afterspikewin',            3 ...
            );
    otherwise
        error('unknown imaging type');
end