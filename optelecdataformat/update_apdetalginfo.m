function ainfo = update_apdetalginfo(ainfo)
apdetalginfo_version = 0.01;
assert(ainfo.version > apdetalginfo_version, 'file version is too new, please update your software');
if ainfo.version < 0.11
    if strcmpi(ainfo.name,'Dynamic programming / nonconvex least squares');
        ainfo.shortname = 'DPNCLS';
    end
end
ainfo = orderfields(ainfo);
ainfo.version = apdetalginfo_version;