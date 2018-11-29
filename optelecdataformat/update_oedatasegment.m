function d = update_oedatasegment(d)
oedatasegment_version = get_oedatasegment_version();
assert(d.version <= oedatasegment_version, 'file version is too new, please update your software');
d.info = update_oedatasegmentinfo(d.info);
if d.version < 0.02
    
    d.imagefile = '';
	d.ephysfile = '';
    
end
if d.version < 0.03
    
    d.frameindex = nan(1, numel(d.t));
    
end
if d.version < 0.04
    
    d.imagefilepartialpath = '';
    d.ephysfilepartialpath = '';    
    
end
if d.version < 0.05
    
    d.f_removed = zeros(size(d.f), class(d.f));
    
end
if d.version < 0.06
    
    d.f_mask = full(d.f_mask);
    
end
if d.version < 0.07
    
    d.roiuid = [];
    
end
if d.version < 0.08
    
    d.f = d.f(:);
    d.f_removed = d.f_removed(:);
    
end
unknown_fields = setdiff(fieldnames(d), fieldnames(empty_oedatasegment));
if ~isempty(unknown_fields)
    
    error(['invalid field: ' unknown_fields{1}]);
    
end
d = orderfields(d);
d.version = oedatasegment_version;