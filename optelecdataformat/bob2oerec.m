function [oerec, info] = bob2oerec(target, opts)
%oerec = bob2oerec(target, opts)

info=struct();
info.roi_idxs = []; % Please implement this!!!

%fixme this function should support cell array of strings

defaultopts = struct( ...
     'minT', 10 ... % ,minimum
     ,'partialpathnames', true ...
     ,'indicator', '' ...
     ,'extractmultiple', false ...
     ,'uncategorized_as_neuron', false ...
     ,'storeimages', false ....
     ,'cellattachedonly', false ...
     ,'mustcontain', '' ...
);
if ~exist('opts','var')
    
    opts=struct;
    
end
opts = mergeOpts(opts,defaultopts, true);
bobloadopts = struct('keepoldalign',false,'minbobversion',1.47);

if isa(target,'struct')
    
    oerec = empty_oerec;
    assert(isfield(target,'image'), 'not a valid bob structure');
    bobdata = target;
    assert(bobdata.version >= bobloadopts.minbobversion); %FIXME update as well?
    if isempty(bobdata.file.bobfile)
        oerec.name = bobdata.image.filename.getFullPath();
    else
        oerec.name = bobdata.file.bobfile;
    end
    oerec.data = bob2oerecdata(bobdata, opts);
    if opts.extractmultiple
        
        %hack, fixme
        n = numel(oerec(1).data);
        oerec = repmat(oerec, 1, n);
        for u = 1:n
            
            oerec(u).data = oerec(u).data(u);
            
        end
        
    end
    
elseif isa(target,'char')
    
    if exist(target, 'file') == 2
        
        bobdata = load_bobdata(target, bobloadopts);
        close_bobdata(target); %don't need to stream from file etc. right now
        oerec = bob2oerec(bobdata, opts);
        return;
        
    elseif exist(target, 'dir')
        assert(~opts.extractmultiple, 'multiple ROIs can be extracted from a bob file only for a single file or struct, not for a directory');
        
        olddir = cd;
        cd(target);
        target = cd; %full pathname, not partial
        cd(olddir);
        bfiles = {}; bdirs = {};
        if target(end) ~= filesep
            target = [target filesep];
        end
        basedir = target;
        nextdir = {target};
        while ~isempty(nextdir)
            
            d = dir(nextdir{1});
            id = cat(2, d.isdir);
            dn = {d.name};
            for u = find(id)
                if strcmpi(dn{u},'.') || strcmpi(dn{u},'..')
                    continue;
                end
                nextdir{end+1} = [nextdir{1} dn{u} filesep]; %#ok<*AGROW>
            end
            b = dir([nextdir{1} '*.bob']);
            for u = 1:numel(b)
                bfiles{end+1} = [nextdir{1} b(u).name];
                bdirs{end+1} = nextdir{1};
            end
            nextdir(1) = [];
            
        end
        oerec = empty_oerec;
        oerec(:) = [];
        [bdirs_unique,~,nind] = unique(bdirs);
                
        for n = 1:max(nind)
            
            if ~isempty(opts.mustcontain) && isempty(strfind(bdirs_unique{n}, opts.mustcontain))
                
                continue;
                
            end
            
            nextneuron_oerec = empty_oerec;
            if opts.partialpathnames
                
                nextneuron_oerec.name = bdirs_unique{n}(numel(basedir) + 1:end);
                
            else
                
                nextneuron_oerec.name = bdirs_unique{n};
                
            end
            if ~isempty(nextneuron_oerec.name) && nextneuron_oerec.name(end) == filesep
                
                nextneuron_oerec.name(end) = [];
                
            end
            for bind = reshape(find(nind == n), 1, [])
                
                try
                    bobdata = load_bobdata(bfiles{bind}, bobloadopts);
                catch ex
                    warning('file %s not imported: %s', bfiles{bind}, ex.message);
                    continue;
                end
                close_bobdata(bobdata); %don't need to stream etc.
                bobdata.file.bobfile = bfiles{bind};
                try
                    
                    nextdata = bob2oerecdata(bobdata, opts);
                    for s = 1:numel(nextdata)
                        
                        if opts.partialpathnames
                            
                            nextdata(s).imagefilepartialpath = bfiles{bind}(numel(basedir) + 1:end);
                            
                        end
                        
                    end
                    nextneuron_oerec.data = [nextneuron_oerec.data nextdata]; %#ok<*NODEF>
                    
                catch ex
                    warning(['Failed to process ' bfiles{bind} ':' ex.message]);
                end
                
            end
            if isempty(nextneuron_oerec.data)
                disp(['No data for neuron in directory ' bdirs_unique{n} ', omitting from results']);
            else
                oerec(1, end + 1) = nextneuron_oerec;
            end
            
        end
        
    else
        
        error('unrecognized input');
        
    end
else
    
    error('unrecognized input');
    
end
if ~isempty(oerec)
    
    oerec = oerec_split_validF(oerec,opts.minT);
    
end
fprintf('\nConverted %d neurons, %d total segments.\n\n',numel(oerec), numel(cat(2, oerec.data)));


function data = bob2oerecdata(bobdata, opts)
data = empty_oedatasegment; data(:) = [];
if isempty(bobdata.image.time) || ~bobdata.roi.numrois
    
    fprintf('no usable data in file %s, skipping\n', bobdata.file.bobfile);
    return;
    
end

if opts.cellattachedonly
    
    carois = reshape(ismember(lower(bobdata.roi.type), {'neuron w/cell-attached'}), 1, []);
    
elseif opts.uncategorized_as_neuron
    
    carois = reshape(ismember(lower(bobdata.roi.type), {'neuron', 'neuron w/cell-attached',''}), 1, []);
    
else
    
    carois = reshape(ismember(lower(bobdata.roi.type), {'neuron', 'neuron w/cell-attached'}), 1, []);
    
end

if bobdata.align.inuse
    
    maplist = setdiff(reshape(unique(bobdata.roi.map(carois)), 1, []), 0);
    
else
    
    maplist = 0;
    
end

bob_roi_each_seg = [];

for m = maplist
    
    if bobdata.align.inuse
        
        fr = find(bobdata.align.map_assoc == m & bobdata.image.goodframes);
        
    else
        
        fr = find(bobdata.image.goodframes);
        
    end
    if isempty(fr), continue; end
    
    maproi = reshape(find(bobdata.roi.map == m & carois), 1, []);
    
    if isempty(maproi)
        
        if bobdata.align.inuse
            
            fprintf('no neuron in map %d of file %s, skipping\n', m, bobdata.file.bobfile);
        else
            
            fprintf('no neuron in file %s, skipping\n', bobdata.file.bobfile);
            
        end
        continue;
        
    end
    
    if numel(maproi) > 1 && ~opts.extractmultiple
        
        if isempty(bobdata.file.bobfile)
            bfs = '';
        else
            bfs = sprintf(' for file %s', bobdata.file.bobfile);
        end
        if bobdata.align.inuse
            warning('%d ROIs on map %d%s, using the first one', numel(maproi), m, bfs);
        else
            warning('%d ROIs%s, using the first one', numel(maproi), bfs);
        end
        maproi = maproi(1);
        
    end
    
    assert(all(maproi <= size(bobdata.roi.kdata, 2)), 'kinetics must be calculated');
    
    for k = maproi
        
        nextdata = bob2oerecdata_singlemap(bobdata, k, opts);
        data = [data nextdata];
        bob_roi_each_seg = [bob_roi_each_seg k];
        
    end
    
end
[~, sortind] = sort(bob_roi_each_seg);
data = data(sortind);


function data = bob2oerecdata_singlemap(bobdata, roiind, opts)
dt = bobdata.image.binsize / 1000;
linet = dt / bobdata.image.height;
pixt = dt / (bobdata.image.width * bobdata.image.height);

centpos = mean(bobdata.roi.interior{roiind},1); %FIXME check basex/basey at centpos if aligned
relt = centpos(1) * pixt + (centpos(2) - 1) * linet - pixt / 2;
t = bobdata.image.time + relt;

data = empty_oedatasegment;
assert(~isempty(bobdata.file.bobfile), 'file must have been saved in .bob format');
data.imagefile = bobdata.file.bobfile;
[~, ff, ee] = fileparts(bobdata.file.bobfile);
data.imagefilepartialpath = [ff ee]; %this may be replaced with a full partial path / file name in the caller

%extract f, performing offset correction if offsets are available
data.f = bobdata.roi.kdata(:, roiind);
darkoffset = bobdata.image.offset;
if isempty(darkoffset) || isnan(darkoffset)
    data.info.darkoffsetcorrected = false;
    data.info.darkoffset = NaN;
else
    data.info.darkoffsetcorrected = true;
    data.info.darkoffset = darkoffset;
    data.f = bobdata.roi.kdata(:, roiind) - darkoffset;
end

data.f_mask = full(~isnan(data.f) & bobdata.image.goodframes & ~bobdata.roi.cellrej(:, roiind));
data.frameindex = 1:bobdata.image.numframes;
data.f_numerictype = 'double';
data.f_offset = 0;
data.f_scale_factor = 1;
data.f_units = 'GSV';
data.t = t;
data.t_offset = 0;
data.t_scale_factor = 1;
data.t_numerictype = 'double';
data.t_units = 'second';

if bobdata.electro.nrOfChannels > 0 % && ~isempty(bobdata.electro.channels{1}.data)
    data.spiketimes = bobdata.electro.channels{1}.data(:,1);
    data.spiketimes_present = true;
    data.spiketimes_window = bobdata.image.time([1 end]) + [1; 1] * dt;
    ns = histc(data.spiketimes, [data.t(1) - dt; data.t]);
    data.apcounts = reshape(ns(1:end-1), numel(data.t), 1);
    data.apcounts_present = true;
end
[~,~,ext] = fileparts(bobdata.image.filename.getFullPath()); %this refers to the imported file, even if a .bob file was later saved
data.info.dwelltimeperimage = dt * size(bobdata.roi.interior{roiind},1) / (bobdata.image.width * bobdata.image.height); %fixme should take into account dead time etc.
if strcmpi(ext,'.cfd')
    
    data.info.acquisitionsoftware = 'CFNT';
    data.info.acquisitionsoftwareversion = str2double(bobdata.image.headerdata.version);
    
    if isfield(bobdata.image.headerdata, 'xpgsettings') && ~isempty(bobdata.image.headerdata.xpgsettings)
        data.info.samplesperpixel = bobdata.image.headerdata.xpgsettings.pct * bobdata.image.headerdata.xpgsettings.cpp * 2;
        data.info.nA2D = size(bobdata.roi.interior{roiind}, 1) * data.info.samplesperpixel;
        
        xpgscalefac = bobdata.image.headerdata.xpgsettings.mult1 * 2 ^ bobdata.image.headerdata.xpgsettings.range1; %multiply this by pixel values to get values in the 20-bit accumulator of 12-bit inputs
        %note that a smaller range/mult give a narrower LUT and larger pixel values for the same inputs
        
        data.info.data_import_scale_factor = xpgscalefac / data.info.samplesperpixel; %this will convert our fluorescence measurements to common units of a2d outputs
        data.f = data.f * data.info.data_import_scale_factor;
    else
        data.info.data_import_scale_factor = NaN;
        warning('xpg settings are not available, using clockfrequency / dwell time to approximate number of A2D samples: %s', data.imagefilepartialpath);
        data.info.nA2D = bobdata.image.headerdata.daClockFrequency * data.info.dwelltimeperimage;
        data.info.samplesperpixel = bobdata.image.headerdata.daClockFrequency * dt / (bobdata.image.width * bobdata.image.height);
    end
    
    data.info.zoomfactor = bobdata.image.headerdata.zoom;
    
    %even though we have a digital offset to compensate the xpg settings in the file, we don't have the full analog / digital offset unless we close the shutter, use a truly dark blood vessel, etc.
    data.info.darkoffsetcorrected = false;
    
elseif isfield(bobdata.image.headerdata,'scanimageheader')
    
    data.info.acquisitionsoftware = 'Scanimage';
    if numel(bobdata.image.headerdata.scanimageheaderstring) > 6 && strcmpi(bobdata.image.headerdata.scanimageheaderstring(1:6),'state.') %3.x
        data.info.acquisitionsoftwareversion = bobdata.image.headerdata.scanimageheader.software.version;
        if bobdata.image.headerdata.scanimageheader.internal.averageSamples
            data.info.data_import_scale_factor = 1;
        else
            data.info.data_import_scale_factor = 1 / bobdata.image.headerdata.scanimageheader.acq.binFactor;
            data.f = data.f / bobdata.image.headerdata.scanimageheader.acq.binFactor; %this puts the data into common units across multiple files that are linearly related to photons / A2D sample (usual caveats apply)
        end
        data.info.darknoise_sd_per_A2D = bobdata.image.headerdata.scanimageheader.acq.pmtOffsetStdDevChannel1; %whether or not to use this value, or how saturation effects ought to be corrected, is a problem for later.
        data.info.samplesperpixel = bobdata.image.headerdata.scanimageheader.acq.binFactor;
        data.info.nA2D = size(bobdata.roi.interior{roiind}, 1) * bobdata.image.headerdata.scanimageheader.acq.binFactor;
        data.info.dwelltimeperimage = data.info.nA2D / bobdata.image.headerdata.scanimageheader.acq.inputRate; %more accurate
        data.info.zoomfactor = bobdata.image.headerdata.scanimageheader.acq.zoomFactor; %should we multiply by basezoomfactor as well?
    elseif isfield(bobdata.image.headerdata, 'scanimageversion') && bobdata.image.headerdata.scanimageversion >= 5 %5.x
        data.info.samplesperpixel = NaN; %no idea whether number of samples per pixel is the same across any given set of files. but we should be able to get it as per vidrio instructions
        data.info.nA2D = size(bobdata.roi.interior{roiind}, 1); %wrong, FIXME
        data.info.data_import_scale_factor = NaN;
        data.info.acquisitionsoftwareversion = bobdata.image.headerdata.scanimageversion;
        if isfield(bobdata.image.headerdata.scanimageheader,'SI5')
            data.info.zoomfactor = bobdata.image.headerdata.scanimageheader.SI5.zoomFactor;
        else
            data.info.zoomfactor = 1;
            warning('Set nonsense value for zoom factor, fix implementation');
        end
    else %4.x
        data.info.acquisitionsoftwareversion = 4; %not clear if this is actually stored in the Thorlabs header or not
        data.info.samplesperpixel = NaN;
        data.info.nA2D = NaN; %no idea whether number of samples per pixel is the same across any given set of files.
        data.info.data_import_scale_factor = NaN;
        %data.info.nA2D = size(bobdata.roi.interior{roiind}, 1); % no idea whether this is correct, probably needs to be adjusted by the number of a2d samples per pixel. in scanimage 4, pixels are sampled at 4-8 MHz depending how you calculate it (i.e. logically or by believing the header)
    end
    
elseif isfield(bobdata.image.headerdata, 'mdfheader')
    
    try %#ok<TRYNC>
        pclock = str2double(bobdata.image.headerdata.mdfheader.Pixel_Clock(1:find(bobdata.image.headerdata.mdfheader.Pixel_Clock == ' ',1) - 1));
        spos = find(bobdata.image.headerdata.mdfheader.Pixel_Clock == 's');
        assert(numel(spos) == 1, 'invalid string');
        switch pcclock(spos - 1)
            case 'p'
                pclock = pcclock * 1e-12;
            case 'n'
                pclock = pcclock * 1e-9;
            case 'u'
                pclock = pcclock * 1e-6;
            case 'm'
                pclock = pcclock * 1e-3;
        end
        data.info.dwelltimeperimage = size(bobdata.roi.interior{roiind}, 1) * pclock;
    end
    data.info.samplesperpixel = NaN;
    data.info.nA2D = size(bobdata.roi.interior{roiind}, 1); %wrong, FIXME
    try %#ok<TRYNC>
        magstring = bobdata.image.headerdata.mdfheader.Magnification;
        assert(~isempty(magstring));
        data.info.zoomfactor = str2double(magstring(1:find(magstring=='x',1)-1));
    end
    
end
data.info.indicator = opts.indicator;
data.info.digitalsettingsindex = 1;
if opts.storeimages
    
    iG = double(bobdata.image.avframe);
    if bobdata.image.nchans > 1
        iR = double(bobdata.image.avframe2);
    else
        iR = zeros(size(iG));
    end
    iB = zeros(size(iG));
    data.image = cat(3, iR, iG, iB);
    
else
    
    data.image = [];
    
end
if bobdata.roi.map(roiind)
    
    data.roi = zeros(bobdata.align.mapdims(bobdata.roi.map(roiind), :));
    
else
    
    data.roi = zeros(bobdata.image.height, bobdata.image.width);
    
end
data.roi(bobdata.roi.interior{roiind}(:,2) + (bobdata.roi.interior{roiind}(:,1) - 1) * size(data.roi, 1)) = 1;
data.roiuid = bobdata.roi.roiIds{roiind};