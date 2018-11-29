function [params,opts] = assign_default_settings_and_params(params,opts,indicatorstring)
nseg = numel(params);
assert(isa(params,'struct') && isa(opts,'struct'), 'invalid inputs');

%first assign indicator-specific defaults where available, then assign remaining fields based on generic defaults:
[dP_indicator, dO_indicator] = sbm.init.params_and_opts_by_indicator(indicatorstring);
[dV, dP] = sbm.init.default_settings_and_params;
fn = union(fieldnames(dV), fieldnames(dO_indicator));

for ii = 1:numel(fn)
    
    if (isfield(dV, fn{ii}) && isa(dV.(fn{ii}), 'struct')) || ...
            (isfield(dO_indicator, fn{ii}) && isa(dO_indicator.(fn{ii}), 'struct'))
        if ~isfield(opts, fn{ii})
            opts.(fn{ii}) = struct();
        end
        subfn = {};
        if isfield(dV, fn{ii})
            subfn = [subfn reshape(fieldnames(dV.(fn{ii})),1,[])]; %#ok<AGROW>
        end
        if isfield(dO_indicator, fn{ii})
            subfn = [subfn reshape(fieldnames(dO_indicator.(fn{ii})),1,[])]; %#ok<AGROW>
        end
        subfn = unique(subfn);
        for jj = 1:numel(subfn)
            if ~isfield(opts.(fn{ii}), subfn{jj})
                if isfield(dO_indicator, fn{ii}) && isfield(dO_indicator.(fn{ii}), subfn{jj})
                    opts.(fn{ii}).(subfn{jj}) = dO_indicator.(fn{ii}).(subfn{jj});
                elseif isfield(dV, fn{ii}) && isfield(dV.(fn{ii}), subfn{jj})
                    opts.(fn{ii}).(subfn{jj}) = dV.(fn{ii}).(subfn{jj});
                end
            end
        end
        opts.(fn{ii}) = orderfields(opts.(fn{ii}));
    else
        if ~isfield(opts, fn{ii})
            if isfield(dO_indicator, fn{ii})
                opts.(fn{ii}) = dO_indicator.(fn{ii});
            else
                opts.(fn{ii}) = dV.(fn{ii});
            end
        end
    end
    
end

fn = union(fieldnames(dP), fieldnames(dP_indicator));
for ii = 1:numel(fn)
    if ~isfield(params, fn{ii})
        for u = 1:nseg
            if isfield(dP_indicator, fn{ii})
                params(u).(fn{ii}) = dP_indicator.(fn{ii});
            else
                params(u).(fn{ii}) = dP.(fn{ii});
            end
        end
    end
end

if ~isfield(opts, 'ResampleThreshold') || isempty(opts.ResampleThreshold)
    opts.ResampleThreshold = opts.Nparticles / 2;
end

if ~isfield(opts, 'Nparticles_prerun') || isempty(opts.Nparticles_prerun)
    opts.Nparticles_prerun = opts.Nparticles;
end
if ~isfield(params, 'fr') && isfield(params, 'gamma_kt_FR')
    for u = 1:nseg
        params(u).fr = prod(params(u).gamma_kt_FR);
    end
end
if ~isfield(params, 'S') && isfield(params, 'normal_meancov_logSR')
    for u = 1:nseg
        params(u).S = exp(params(u).normal_meancov_logSR(1, 1));
    end
end
if ~isfield(params, 'R') && isfield(params, 'normal_meancov_logSR')
    for u = 1:nseg
        params(u).R = exp(params(u).normal_meancov_logSR(2, 1));
    end
end

params = orderfields(params);
opts = orderfields(opts);