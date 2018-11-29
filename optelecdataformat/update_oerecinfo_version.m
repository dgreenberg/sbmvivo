function oerecinfo = update_oerecinfo_version(oerecinfo)
oerecinfo_version = 0.012;
assert(oerecinfo.version <= oerecinfo_version, 'file version is too new, please update your software');
if oerecinfo.version < 0.011
    oerecinfo.attributes_qual = {};
    oerecinfo.attributes_quant = {};
    oerecinfo.values_quant = []; %double precision
end
if oerecinfo.version < 0.012
    if isfield(oerecinfo,'attributes_quan')
        oerecinfo.attributes_quant = oerecinfo.attributes_quan;
        oerecinfo = rmfield(oerecinfo, 'attributes_quan');
    end
    if iscell(oerecinfo.attributes_quant) && numel(oerecinfo.attributes_quant) && iscell(oerecinfo.attributes_quant{1})
        oerecinfo.attributes_quant = {};
    end
    if iscell(oerecinfo.attributes_qual) && numel(oerecinfo.attributes_qual) && iscell(oerecinfo.attributes_qual{1})
        oerecinfo.attributes_qual = {};
    end
end
unknown_fields = setdiff(fieldnames(oerecinfo), fieldnames(empty_oerecinfo));
if ~isempty(unknown_fields)
    error(['invalid field: ' unknown_fields{1}]);
end
oerecinfo = orderfields(oerecinfo);
oerecinfo.version = oerecinfo_version;