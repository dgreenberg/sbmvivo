function [st, w, twin] = convert_spikecounts_to_weighted_spiketrain(ns, t_ns)
assert(~all(isnan(ns)), 'spikecounts cannot all be NaN');

s = find(~isnan(ns), 1);
e = find(~isnan(ns), 1, 'last');
ns = ns(s:e);
t_ns = t_ns(s:e);

[index, ~, w] = find(ns);
st = t_ns(index);
dt = median(diff(t_ns));
twin = [t_ns(1) - dt / 2, t_ns(end) + dt / 2];