function simvdata = sbm_simavspikeresponse_eachneuron(oedb, fileinfo, ds, algind, avwin, dt_sim, T_sim_sec)
if ~exist('avwin', 'var') || isempty(avwin)
    
    avwin = [0 1];
    
end
if ~exist('dt_sim', 'var') || isempty(dt_sim)
    
    dt_sim = 1e-3;

end
if ~exist('T_sim_sec', 'var') || isempty(T_sim_sec)

    T_sim_sec = 3;    
    
end
mintbefore_lowfr = 1.5; %sec
mintbefore_highfr = 0.3; %sec
frthresh = 1; %Hz
T_sim = ceil(T_sim_sec / dt_sim);
oerec = fetch_dataset_oerecarray(oedb, fileinfo, ds);
nspikes = 1;

totalT = arrayfun(@(x) sum(arrayfun(@(d) max(d.t) - min(d.t), x.data)), oerec);
totalaps = arrayfun(@(x) sum(arrayfun(@(d) sum(d.apcounts), x.data)), oerec);
fr = totalaps ./ totalT;

mintbefore = mintbefore_lowfr * ones(1, numel(oerec));
mintbefore(fr > frthresh) = mintbefore_highfr;
[avval_data, avval_sim, peakval_data, peakval_sim, nsingles] = deal(nan(1, oedb.nneurons(ds)));
[dff_sim, dff_data, t_sim, t_data] = deal(cell(1, oedb.nneurons(ds)));

for n = 1:oedb.nneurons(ds)
    
    dt_eachseg = arrayfun(@(v) diff(v.t(1:2)) * v.t_scale_factor, oerec(n).data);
    
    maxdt_thisneuron = max(dt_eachseg);
    
    results = fetch_results(oedb, fileinfo, ds, n, 1, algind);  % first segment    
    if isempty(fieldnames(results.params)), continue; end
    
    [params, sbm_opts] = sbm.init.assign_default_settings_and_params(results.params, results.opts, oerec(n).data(1).info.indicator);    
    
    [brightness, states, t_sim{n}, spike_time] = sbm.model.spikeresponse(params, sbm_opts, dt_sim, T_sim, nspikes);
    t_sim{n} = t_sim{n} - spike_time;
    dff_sim{n} = (brightness - brightness(1)) / brightness(1);
    
    [av, t_data{n}, ~, ~, ~, nap_isolated] = avspikeresponse(oerec(n), mintbefore(n), maxdt_thisneuron, true);
    dff_data{n} = av{1};
    nsingles(n) = nap_isolated(1);
    
    avval_data(n) = mean(av{1}(t_data{n} > avwin(1) & t_data{n} < avwin(2)));
    avval_sim(n) =  mean(dff_sim{n}(t_sim{n} > avwin(1) & t_sim{n} < avwin(2)));
    
    peakval_data(n) = max(av{1}(t_data{n} > avwin(1) & t_data{n} < avwin(2)));
    peakval_sim(n) =  max(dff_sim{n}(t_sim{n} > avwin(1) & t_sim{n} < avwin(2)));
    
end
simvdata = orderfields(struct( ...
    'avwin', avwin, ...
    'dt_sim', dt_sim, ...
    'T_sim_sec', T_sim_sec, ...
    'T_sim', T_sim, ...
    'dff_data', {dff_data}, ...
    'dff_sim', {dff_sim}, ...
    't_data', {t_data}, ...
    't_sim', {t_sim}, ...
    'avval_data', avval_data, ...
    'avval_sim', avval_sim, ...
    'peakval_data', peakval_data, ...
    'peakval_sim', peakval_sim, ...
    'nsingles', nsingles ...
    ));
