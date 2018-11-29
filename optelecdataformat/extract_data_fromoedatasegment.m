function [f,it,dt,ns,indicatorstring,nA2D,st,f_removed,stwin] = extract_data_fromoedatasegment(d)
%[f,it,dt,ns,indicatorstring,nA2D,st,f_removed] = extract_data_fromoedatasegment(d)
f = reshape(double(d.f) * d.f_scale_factor + d.f_offset,[],1);
f_removed = double(d.f_removed) * d.f_scale_factor + d.f_offset;
if ~isempty(f_removed)
    
    assert(size(f_removed, 1) == numel(f), 'invalid shape for f_removed');
    
end
it = reshape(double(d.t) * d.t_scale_factor + d.t_offset,[],1);
dt = median(diff(it));

[ns, st] = deal([]);
stwin = nan(1, 2);

if d.apcounts_present

    ns = d.apcounts;
    
end

if d.spiketimes_present

    st = d.spiketimes; %filter by window? FIXME
    stwin = d.spiketimes_window(:)';
    
end

switch  d.f_units
    
    case 'df/f0'
        
        f = f * 100 + 100;
        
    case 'GSV'
        
        %do nothing
        
    otherwise
        
        warning('oedatabrowser:unitsnotrecognized','unrecognized fluorescence units');
        
end

indicatorstring = lower(d.info.indicator);
nA2D = d.info.nA2D;

