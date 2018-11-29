function [settingsindex, sessionindex, darkoffsetcorrected, fdcindex, focalplaneindex] = oerec_datagroupindices_singleneuron(oerec)
assert(numel(oerec) == 1, 'one neuron at a time');
if any(cellfun(@isempty, {oerec.data.imagefilepartialpath}))
    warning('some partial paths are empty, cannot determine whether data is from the same file or not. all data without partial paths will be treated as if from the same file');
end
[~, ~, fileindex] = unique({oerec.data.imagefilepartialpath});
allinfo = cat(2, oerec.data.info);
settingsindex = cat(1, allinfo.settingsindex);
focalplaneindex = cat(1, allinfo.focalplaneindex);
darkoffsetcorrected = cat(1, allinfo.darkoffsetcorrected);
[~, ~, fdcindex] = unique(settingsindex);
fdcindex(darkoffsetcorrected) = 0;
[~, ~, sessionindex] = unique([settingsindex(:) focalplaneindex(:) fileindex(:)],'rows'); %at present always use different F0 for different files or for parts of the same file with different fpis (e.g. multiple motion-correction maps)