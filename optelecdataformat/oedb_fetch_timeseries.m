function [f,it,dt,ns,indicatorstring,nA2D,st,f_removed,oerec,stwin] = oedb_fetch_timeseries(oedb, fileinfo, dataset, neuron, segment) %f,it column vectors; ns row vector
oerec = fetch_neurondata(oedb, fileinfo, dataset, neuron);
[f,it,dt,ns,indicatorstring,nA2D,st,f_removed,stwin] = extract_data_fromoerec(oerec, segment);
f = f{1}; it = it{1}; dt = dt(1); ns = ns{1}; indicatorstring = indicatorstring{1}; nA2D = nA2D(1); f_removed = f_removed{1}; stwin = stwin(1, :);
if isempty(st)
    st = [];
else
    st = st{1};
end