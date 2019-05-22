function nname = oedb_generate_neuron_name(oerec, v)
if isempty(oerec.name)
    nname = ['nrn' num2str(v)];
else
    nname = oerec.name;
end
attr = reshape(oerec.info.attributes_qual, 1, []);
if ~isempty(oerec.data) && ~isempty(oerec.data(1).info.indicator)
    attr = [{oerec.data(1).info.indicator} attr];
end
if ~isempty(oerec.info.species)
    attr = [{oerec.info.species} attr];
end
if numel(oerec.data) > 1
    attr = [attr {[num2str(numel(oerec.data)) ' segments']}];
end
if numel(attr) > 0
    nname = [nname ' (' attr{1}];
    for jj = 2:numel(attr)
        nname = [nname ', ' attr{jj}];
    end
    nname = [nname ')'];
end