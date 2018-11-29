function oedsinfo = update_oedatasegmentinfo(oedsinfo)
oedatasegmentinfo_version = 0.19;
assert(oedsinfo.version <= oedatasegmentinfo_version, 'file version is too new, please update your software');
if oedsinfo.version < 0.11
    oedsinfo.digitalsettingsindex = NaN;
end
if oedsinfo.version < 0.12
    oedsinfo.dwelltimeperimage = NaN;
    oedsinfo.nA2D = NaN;
end
if oedsinfo.version < 0.13
    oedsinfo.settingsindex = NaN;
end
if oedsinfo.version < 0.14
    if isa(oedsinfo.zoomfactor,'char')
        oedsinfo.zoomfactor = NaN;
    end
end
if oedsinfo.version < 0.15
    oedsinfo.focalplaneindex = NaN;
end
if oedsinfo.version < 0.16
    oedsinfo.darkoffsetcorrected = false;
    oedsinfo.darkoffset = NaN;
end
if oedsinfo.version < 0.17
    oedsinfo.data_import_scale_factor = NaN;
end
if oedsinfo.version < 0.18
    oedsinfo.darknoise_sd_per_A2D = NaN;
end
if oedsinfo.version < 0.19
    oedsinfo.samplesperpixel = NaN;
end
unknown_fields = setdiff(fieldnames(oedsinfo), fieldnames(empty_oedatasegmentinfo));
if ~isempty(unknown_fields)
    error(['invalid field: ' unknown_fields{1}]);
end
oedsinfo = orderfields(oedsinfo);
oedsinfo.version = oedatasegmentinfo_version;