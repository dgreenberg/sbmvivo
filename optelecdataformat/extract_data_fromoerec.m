function [f,it,dt,ns,indicatorstring,nA2D,st,f_removed,stwin] = extract_data_fromoerec(oerec, segmentlist)
assert(numel(oerec) == 1,'this function is for one neuron at a time');
if ~exist('segmentlist','var')
    
    segmentlist = 1:numel(oerec.data);
    
else
    
    assert(~isempty(segmentlist), 'segmentlist cannot be empty');
    
end
[f, f_removed, it, ns, indicatorstring,st] = deal(cell(1, numel(segmentlist)));
[dt, nA2D] = deal(nan(1, numel(segmentlist)));
stwin = nan(numel(segmentlist), 2);
for k = 1:numel(segmentlist)
    
    segment = segmentlist(k);    
    [f{k}, it{k}, dt(k), ns{k}, indicatorstring{k}, nA2D(k), st{k}, f_removed{k}, stwin(k, :)] = ...
        extract_data_fromoedatasegment(oerec.data(segment));
    
end