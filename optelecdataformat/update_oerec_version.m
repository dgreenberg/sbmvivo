function oerec = update_oerec_version(oerec)
if numel(oerec) > 1
    newoerec = repmat(empty_oerec, 1, numel(oerec));    
    for u = 1:numel(oerec)
        newoerec(u) = update_oerec_version(oerec(u));
    end
    oerec = newoerec;
    return;
end
oerec_version = 0.01;
assert(oerec.version <= oerec_version, 'file version is too new, please update your software');
oerec.info = update_oerecinfo_version(oerec.info);
data_new = repmat(empty_oedatasegment, 1, numel(oerec.data));
for n = 1:numel(oerec.data)
    data_new(n) = update_oedatasegment(oerec.data(n));
end
oerec.data = data_new;
%if oerec.version < XXX
    %update fields
%end
unknown_fields = setdiff(fieldnames(oerec), fieldnames(empty_oerec));
if ~isempty(unknown_fields)
    error(['invalid field: ' unknown_fields{1}]);
end
oerec = orderfields(oerec);
oerec.version = oerec_version;