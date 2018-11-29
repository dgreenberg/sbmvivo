function oesim = update_oesim_version(oesim)
if numel(oesim) > 1
    newoesim = repmat(empty_oesimulation, 1, numel(oesim));
    for u = 1:numel(oesim)
        newoesim(u) = update_oesim_version(oesim(u));
    end
    oesim = newoesim;
    return;
end
oesim_version = 0.01;
assert(oesim.version <= oesim_version, 'file version is too new, please update your software');
oesim.oerec = update_oerec_version(oesim.oerec);
%if oesim.version < XXX
    %update fields
%end
oesim = orderfields(oesim);
oesim.version = oesim_version;