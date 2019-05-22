function s = samedatasizemat(oedb)
s = false(oedb.ndatasets);
for u1 = 1:oedb.ndatasets
    for u2 = u1 + 1:oedb.ndatasets
        if oedb.nneurons(u1) ~= oedb.nneurons(u2), continue; end
        if any(oedb.nsegments{u1} ~= oedb.nsegments{u2}), continue; end
        d1 = cat(2, oedb.data{u1}.data);
        d2 = cat(2, oedb.data{u2}.data);
        if any(cellfun(@numel, {d1.f}) ~= cellfun(@numel, {d2.f})), continue; end
        s(u1, u2) = true; s(u2, u1) = true;
    end
end