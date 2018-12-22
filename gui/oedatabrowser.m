function varargout = oedatabrowser(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @oedatabrowser_OpeningFcn, ...
    'gui_OutputFcn',  @oedatabrowser_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


function oedatabrowser_OpeningFcn(hObject, eventdata, handles, varargin)
if strcmpi(get(hObject,'visible'),'on') %already running
    figure(hObject);
    return;
end
handles.output = hObject;
guidata(hObject, handles);
%should read an ini file here FIXME
fileinfo = default_oedatabrowser_fileinfostruct;
setappdata(handles.oedatabrowser_figure,'fileinfo',fileinfo);

try
    oedb_cleartmpdir(fileinfo);
catch ex
    uiwait(errordlg(ex.message, 'Cannot clear tmp directory', 'modal'));
    delete(handles.oedatabrowser_figure);
    return;
end

try
    oedb = oedb_init(fileinfo);
catch ex
    uiwait(errordlg(ex.message, 'Initialization failed', 'modal'));
    delete(handles.oedatabrowser_figure);
    return;
end

if ~oedb.opt.storedatainmemory
    oedb = oedb_movedatatodisk(oedb, fileinfo);  %fixme is this needed?
end

setappdata(handles.oedatabrowser_figure,'oedb',oedb);
setappdata(handles.oedatabrowser_figure,'ui',default_oedatabrowser_uistruct);
update_ui(handles, oedb, fileinfo);
show_data(handles);
init_appearance(handles);
set(hObject,'renderer','painters');  % zbuffer causes hangs when zooming in / out
set(handles.menu_train,'visible','off'); % feature not yet fully implemented


function init_appearance(handles)
vlevel = 3; nvlevels = 10;
set([handles.f_axes handles.spike_axes handles.ca_axes],'tickdir','out');
set([handles.f_axes handles.spike_axes],'xticklabel',[]);
vmenu = nan(1, nvlevels);
for v = 1:nvlevels
    vmenu(v) = uimenu('parent',handles.sbm_verbosity,'label',num2str(v),'callback',@verbositylevel_Callback,'checked','off');
end
set(vmenu(vlevel),'checked','on');


function oedb = oedb_init(fileinfo)
oedb = empty_oedbstruct;
oedb = update_data(oedb, fileinfo);
oedb = update_algs(oedb);
oedb = oedb_init_outputs(oedb);
oedb = oedb_init_stats(oedb);


function oedb = update_algs(oedb)
[oedb.algfunctions, oedb.alginfo, ~, ~, oedb.trainfunctions] = oerec_getalginfo;
oedb.algnames = cell(1, numel(oedb.alginfo));
for k = 1:numel(oedb.alginfo)
    oedb.algnames{k} = oedb.alginfo(k).shortname;
end
[oedb.algavailable] = deal(true(1, numel(oedb.alginfo)));


function oedb = import_dataset(oedb, fileinfo, datafile)
assert(exist(datafile, 'file') == 2, 'file not found');

u = oedb.ndatasets + 1;
oedb.sourcedatafiles{u} = datafile;

nneuronstotal = sum(oedb.nneurons);
nsegtotal = sum(cat(2, oedb.nsegments{:}));

[~, ff, ee] = fileparts(datafile);
if numel(ff) > 0 && strcmpi(ee, '.mat')
    oedb.datasetnames{u} = ff;
else
    oedb.datasetnames{u} = [ff ee];
end

if ~isempty(whos('-file', datafile, 'oerec'))
    
    oedb.issimulation(u) = false;
    if numel(oedb.datasetnames{u}) > 8 && strcmpi(oedb.datasetnames{u}(1:8),'dataset_')
        oedb.datasetnames{u} = oedb.datasetnames{u}(9:end);
    end
    oedb.data{u} = oerec_load(oedb.sourcedatafiles{u}); %also updates version if necessary
    oedb.nneurons(u) = numel(oedb.data{u});
        
elseif ~isempty(whos('-file', datafile, 'oesim'))
    
    oedb.issimulation(u) = true;
    if numel(oedb.datasetnames{u}) > 11 && strcmpi(oedb.datasetnames{u}(1:11),'simulation_')
        oedb.datasetnames{u} = oedb.datasetnames{u}(12:end);
    end
    oedb.data{u} = oesim_load(oedb.sourcedatafiles{u}); %also updates version if necessary
    oedb.nneurons(u) = numel(oedb.data{u}.oerec);
    
else
    
    error('invalid file');
    
end
assert(~isempty(oedb.data{u}), ['failed to load data from file ' oedb.sourcedatafiles{u}]);

oedb.expandneuron{u} = true(1, oedb.nneurons(u));
[oedb.neuronnames{u}, oedb.segmentnames{u}, oedb.segcode{u}, oedb.neuroncode{u}] = deal(cell(1, oedb.nneurons(u)));
oedb.nsegments{u} = nan(1, oedb.nneurons(u));
for v = 1:oedb.nneurons(u)
    
    nneuronstotal = nneuronstotal + 1;
    
    if oedb.issimulation(u)
        oerec = oedb.data{u}.oerec(v);
        oesim = oedb.data{u}; oesim.oerec = oesim.oerec(v);
    else
        oerec = oedb.data{u}(v);
    end
    oedb.nsegments{u}(v) = numel(oerec.data);
    if oedb.nsegments{u}(v) < 2
        oedb.expandneuron{u}(v) = false;
    end
    
    oedb.neuronnames{u}{v} = oedb_generate_neuron_name(oerec, v);
    oedb.neuroncode{u}{v} = num2str(nneuronstotal, 9);
    if ~oedb.opt.storedatainmemory
        if oedb.issimulation(u)
            save([fileinfo.tmpdir 's' oedb.neuroncode{u}{v} '.mat'],'oesim');
        else
            save([fileinfo.tmpdir 'n' oedb.neuroncode{u}{v} '.mat'],'oerec');
        end
    end
    
    ifiles = {oerec.data.imagefilepartialpath};
    ifiles_fullpath = {oerec.data.imagefilepartialpath};
    ifiles(cellfun(@isempty, ifiles)) = ifiles_fullpath(cellfun(@isempty, ifiles));
    for w = 1:oedb.nsegments{u}(v)
        if ~isempty(ifiles{w})
            [~,ifiles{w},~] = fileparts(ifiles{w});
        end
    end
    [~,~,ifileind] = unique(ifiles);
    
    [oedb.segcode{u}{v}, oedb.segmentnames{u}{v}] = deal(cell(1, oedb.nsegments{u}(v)));
    for w = 1:oedb.nsegments{u}(v)
        nsegtotal = nsegtotal + 1;
        oedb.segcode{u}{v}{w} = num2str(nsegtotal, 9);
        if ~isempty(oerec.data(w).imagefile)
            oedb.segmentnames{u}{v}{w} = ifiles{w};
            if sum(ifileind == ifileind(w)) > 1
                ifilesegind = find(w == find(ifileind == ifileind(w)));
                oedb.segmentnames{u}{v}{w} = [oedb.segmentnames{u}{v}{w} ' (seg ' num2str(ifilesegind) ')'];
            end
        else
            oedb.segmentnames{u}{v}{w} = num2str(w);
        end
        nmatchingchars = strnmatchingchars(oedb.segmentnames{u}{v}{w}, oedb.neuronnames{u}{v});
        if nmatchingchars > 5
            oedb.segmentnames{u}{v}{w} = oedb.segmentnames{u}{v}{w}(nmatchingchars + 1:end); %remove neuron name from segment name
        end
        if numel(oedb.segmentnames{u}{v}{w}) > 1 && ismember(oedb.segmentnames{u}{v}{w}(1), '/\_:., ')
            oedb.segmentnames{u}{v}{w} = oedb.segmentnames{u}{v}{w}(2:end);
        end
    end
    assert(nsegtotal < oedb_maxtotalsegs, ['oedatabrowser is restricted to less than ' num2str(oedb_maxtotalsegs) ' data segments']);
    
end
if ~oedb.opt.storedatainmemory
    
    oedb.data{u} = [];
    
end
oedb.ndatasets = u;
oedb.expanddataset(u) = true;
oedb.samedatasize = samedatasizemat(oedb);  % keep track of whether each pair of datasets can be compared fluorscence value by fluorescence value


function oedb = update_data(oedb, fileinfo)
[oerec_datasets, oesim_datasets, oe_datadir] = oerec_getdatainfo();
for u = 1:numel(oerec_datasets)
    
    oedb = import_dataset(oedb, fileinfo, [oe_datadir oerec_datasets{u}]);
    
end
for u = 1:numel(oesim_datasets)
    
    oedb = import_dataset(oedb, fileinfo, [oe_datadir oesim_datasets{u}]);
    
end


function file_import_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles);
[f, p] = getfile_fromp(fileinfo.odbdir, '*.mat', 'Import optical / electrical dataset');
oedb = import_dataset(oedb, fileinfo, [p f]);
oedb = oedb_init_outputs(oedb, oedb.ndatasets);
oedb = oedb_init_stats(oedb, oedb.ndatasets);
setappdata(handles.oedatabrowser_figure,'oedb',oedb);
update_ui(handles, oedb, fileinfo);


function oedb = empty_oedbstruct
oedb = struct('algnames',{{}},'algfunctions',{{}},'trainfunctions',{{}},'alginfo',empty_apdetalginfo,'algavailable',false(1,0), ...
    'sourcedatafiles',{{}},'datasetnames',{{}},'neuronnames',{{}},'segmentnames',{{}},...
    'ndatasets',0,'nneurons',zeros(1,0),'nsegments',{{}},'issimulation',false(1,0),...    
    'expanddataset',false(1,0),'expandneuron',{{}},'opt',oedb_defaultoptstruct,...
    'data',{{}},'results',{{}},'segcode',{{}} ...
    ,'stats', orderfields(struct('bydataset', repmat(empty_statstruct,0,0), 'byneuron', {{}}, 'bysegment', {{}})) ...
    ,'modified',false,'version',oedb_version);
oedb = orderfields(oedb);


function m = oedb_maxtotalsegs
m = 1e9 - 1;


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


function n = strnmatchingchars(s1, s2)
m = min(numel(s1), numel(s2));
n = sum(cumprod(double(s1(1:m) == s2(1:m))));


function ui = default_oedatabrowser_uistruct
ui = struct('row2dataset',[],'row2neuron',[],'row2segment',[],'selecteddataset',0,'selectedneuron',0,'selectedsegment',0);


function zoomoutfull(handles)
[ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles);
[~,it,~,~] = oedb_fetch_timeseries(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment);
dt = median(diff(it));
XL = [it(1) it(end)] + [-2 2] * dt;
set([handles.f_axes handles.spike_axes handles.ca_axes],'xlim',XL);


function display_datainfotext(oerec, segment, axh)
%display some info about this data
zoomfac = oerec.data(segment).info.zoomfactor;
if ~isnan(zoomfac)
    
    zoomstr = ['zoom ' num2str(zoomfac)];
    
else
    
    zoomstr = '';
    
end
acqsoftwarename = oerec.data(segment).info.acquisitionsoftware;
acqsoftwareversion = oerec.data(segment).info.acquisitionsoftwareversion;

[f,it,~,ns,~,~,~,~] = extract_data_fromoerec(oerec, segment);
f = f{1}; it = it{1}; ns = ns{1};

dt = median(diff(it));

texty = (max(f) + min(f)) / 2;

XL = [it(1) it(end)] + [-2 2] * dt;

truefr = mean(ns(~isnan(ns))) / dt;

str = {};
if ~isnan(truefr)
    
    str{end + 1, 1} = sprintf('%g spikes/sec', truefr);
    
end

str{end + 1, 1} = sprintf('%s %g @ %g Hz, %d frames', acqsoftwarename, acqsoftwareversion, 1 / dt, numel(f));
str{end + 1, 1} = zoomstr;

abs_diff_F_median = median(abs(diff(f)));
halfnormal_median = 2 * erfinv(0.5); %median of the absolute difference of two independent normal variables with zero mean and unit variance, see https://en.wikipedia.org/wiki/Half-normal_distribution
total_F_sd = abs_diff_F_median / halfnormal_median; %a median-based estimate of fluorescence noise s.d.
total_F_sd = max(total_F_sd, std(f) / 20);
str{end + 1} = sprintf('total s.d.: %f', total_F_sd);
str{end + 1} = sprintf('zeta est: %g', total_F_sd ^ 2 * oerec.data(segment).info.nA2D);

text(XL(1) + diff(XL) / 30, texty, str,'parent', axh, 'verticalalignment', 'middle');


function display_statsinfotext(oedb, fileinfo, dataset, neuron, segment, algind, selectiontype, axh, it)
dt = median(diff(it));
rheight = 0.5; %max([max(estns) max(ns) 0.5]) * 0.8;

switch lower(selectiontype)
    case 'dataset'
        
        s = oedb.stats.bydataset(algind, dataset);
        
    case 'neuron'
        
        s = oedb.stats.byneuron{dataset}(algind, neuron);
        
    case 'segment'
        
        s = oedb.stats.bysegment{dataset}{neuron}(algind, segment);
        
    otherwise
        
        error('unrecognized selectiontype');
        
end

t = cell(0, 1); %cell array of strings that will displayed as text
fn = {'corr_spikecounts', 'det_spikecounts', 'fp_spikecounts', 'corr_spiketimes', 'det_spiketimes', 'fp_spiketimes'};

for u = 1:numel(fn)
    
    sdval = s.([fn{u} '_sd']);
    
    if isempty(s.(fn{u})) || isnan(s.(fn{u}))
        
        continue;
        
    elseif isempty(sdval) || isnan(sdval)
        
        t{end + 1, 1} = sprintf('%s: %g', fn{u}, s.(fn{u}));
        
    else
        
        t{end + 1, 1} = sprintf('%s: %g +/- %g', fn{u}, s.(fn{u}), s.([fn{u} '_sd']));
        
    end
    
end

XL = [it(1) it(end)] + [-2 2] * dt;
text(XL(1) + diff(XL) / 30, rheight, t, 'parent', axh, 'interpreter', 'none', 'verticalalignment', 'bottom');


function y = estspiketrainy0
y = 1;


function exportplot_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, algind, algname, oedb] = oedb_getbasics(handles);
if ~oedb.ndatasets, return; end
results = fetch_results(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment, algind);
if any([ui.selecteddataset ui.selectedneuron ui.selectedsegment] == 0), return; end;
[f,it,~,ns,~,~,st,f_removed] = oedb_fetch_timeseries(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment); %f,it column vectors; ns row vector
dt = median(diff(it));

figh = figure('renderer', 'painters', 'name', sprintf('%s segment %d of %d', ...
    oedb.neuronnames{ui.selecteddataset}{ui.selectedneuron}, ui.selectedsegment, oedb.nsegments{ui.selecteddataset}(ui.selectedneuron)));
showstates = strcmpi(algname, 'sbm') && isfield(results.opts, 'model');
if showstates
    
    if strcmpi(results.opts.model, '5s2b')
        
        nstates = 6 + numel(results.params.Btot);
        
    elseif strcmpi(opts.model, '3s')
        
        nstates = 3;
        
    end
    f_axes = subplot(nstates + 2, 1, 1, 'parent', figh);
    spike_axes = subplot(nstates + 2, 1, 2, 'parent', figh);
    for k = 1:nstates
        
        state_axes(k) = subplot(nstates + 2, 1, 2 + k);
        
    end
    
else
    
    f_axes = subplot(2, 1, 1, 'parent', figh);
    spike_axes = subplot(2, 1, 2, 'parent', figh);
    state_axes = [];
    
end

plot_f_and_spikes(it, f, st, ns, [f_axes spike_axes], f_removed, ui_ison(handles.showremovedf));
plot_apdetresults(results, dt, f_axes, spike_axes, algname, state_axes, false);
linkaxes([f_axes spike_axes state_axes], 'x');
set([f_axes spike_axes state_axes], 'xlim', [it(1) - dt, it(end) + dt]);


function show_data(handles)
[ui, fileinfo, selectiontype, algind, algname, oedb] = oedb_getbasics(handles);
if ~oedb.ndatasets, return; end
clear_axes_oedatabrowser(handles);
if any([ui.selecteddataset ui.selectedneuron ui.selectedsegment] == 0), return; end;
[f,it,~,ns,~,~,st,f_removed,oerec] = oedb_fetch_timeseries(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment); %f,it column vectors; ns row vector
dt = median(diff(it));
st = st(st >= it(1) - dt & st <= it(end) + dt);

results = fetch_results(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment, algind);

plot_f_and_spikes(it, f, st, ns, [handles.f_axes handles.spike_axes], f_removed, ui_ison(handles.showremovedf));
plot_apdetresults(results, dt, handles.f_axes, handles.spike_axes, algname, handles.ca_axes, true);

display_datainfotext(oerec, ui.selectedsegment, handles.f_axes);
display_statsinfotext(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment, algind, selectiontype, handles.spike_axes, it);


function [model_params, model_data, siminfo] = fetch_siminfo(oedb, fileinfo, dataset, neuron)
assert(oedb.issimulation(dataset), 'dataset must be a simulation');
oesim = fetch_simstruct(oedb, fileinfo, dataset, neuron);
model_params = oesim.model_parameters;
model_data = oesim.model_data;
siminfo = rmfield(oesim,{'model_data','model_parameters','oerec'});

function oesim = fetch_dataset_oesimarray(oedb, fileinfo, dataset)
if oedb.opt.storedatainmemory
    
    oesim = oedb.data{dataset};
    
else
    for v = 1:oedb.nneurons(dataset)
        
        oesim(v) = oesim_load([fileinfo.tmpdir 's' oedb.neuroncode{dataset}{neuron} '.mat']); %also updates version if necessary
        
    end
end


function oesim = fetch_simstruct(oedb, fileinfo, dataset, neuron)
if oedb.opt.storedatainmemory
    
    oesim = oedb.data{dataset}(neuron);
    
else
    
    oesim = oesim_load([fileinfo.tmpdir 's' oedb.neuroncode{dataset}{neuron} '.mat']); %also updates version if necessary
    
end


function oedb = modify_neurondata(oedb, fileinfo, dataset, neuron, oerec)
assert(numel(oerec) == 1, 'this function is for one neuron at a time');
if oedb.issimulation(dataset)
    oesim = fetch_simstruct(oedb, fileinfo, dataset, neuron);
    oesim.oerec = oerec;
    if oedb.opt.storedatainmemory
        
        oedb.data{dataset}(neuron) = oesim;
        
    else
        
        save([fileinfo.tmpdir 's' oedb.neuroncode{dataset}{neuron} '.mat'],'oesim');
        
    end
else
    if oedb.opt.storedatainmemory
        
        oedb.data{dataset}(neuron) = oerec;
        
    else
        
        save([fileinfo.tmpdir 'n' oedb.neuroncode{dataset}{neuron} '.mat'],'oerec');
        
    end
end


function clear_axes_oedatabrowser(handles)
for axh = [handles.f_axes handles.spike_axes handles.ca_axes]
    delete(get(axh,'children'));
end


function opt = oedb_defaultoptstruct %FIXME: should probably turn off storedatainmemory automatically when data on disk corresponds to at least some minimum predicted size in memory
opt = orderfields(struct(...
    'storedatainmemory',true,'figsavedir',[fileparts(mfilename('fullpath')) filesep '..' filesep 'figs' filesep] ...
    ,'expandonlyselected',true ...
    ,'tol_sec', 0.2 ...
    ));


function update_ui(handles, oedb, fileinfo)
update_datalist(handles, oedb);
update_alglist(handles, oedb);
update_ui_oedbopts(handles, oedb);
update_figurename(handles.oedatabrowser_figure, oedb, fileinfo);
update_selection(handles);


function update_figurename(figh, oedb, fileinfo)
if isempty(fileinfo.odbfile)
    
    set(figh,'name','oedatabrowser');
    
elseif oedb.modified
    
    set(figh,'name',['oedatabrowser -- ' fileinfo.odbfile '*']);
    
else
    
    set(figh,'name',['oedatabrowser -- ' fileinfo.odbfile]);
    
end


function update_ui_oedbopts(handles, oedb)
if oedb.opt.storedatainmemory
    
    set(handles.opt_storedatainmemory,'checked','on');
    
else
    
    set(handles.opt_storedatainmemory,'checked','off');
    
end
ui_setval(handles.tol_edit, oedb.opt.tol_sec * 1000);


function update_datalist(handles, oedb, selectiontype)
if ~exist('selectiontype','var')
    selectiontype = get_selectiontype(handles);
end
tabstr = '    ';
ui = getappdata(handles.oedatabrowser_figure,'ui');
listnames = {}; ui.row2dataset = []; ui.row2neuron = []; ui.row2segment = [];
if ui.selecteddataset == 0 && oedb.ndatasets > 0, ui.selecteddataset = 1; end
for u = 1:oedb.ndatasets
    listnames{1, end+1} = [oedb.datasetnames{u} ' (' num2str(oedb.nneurons(u)) ' neurons)']; %#ok<*AGROW>
    if oedb.issimulation(u)
        listnames{1, end} = [listnames{1, end} ' (simulation)']; %#ok<*AGROW>
    end
    ui.row2dataset(1, end+1) = u; ui.row2neuron(1, end+1) = 0; ui.row2segment(1, end+1) = 0;
    if ~oedb.expanddataset(u) || (oedb.opt.expandonlyselected && u ~= ui.selecteddataset), continue; end;
    for v = 1:oedb.nneurons(u)
        listnames{1, end+1} = [tabstr oedb.neuronnames{u}{v}];
        ui.row2dataset(1, end+1) = u; ui.row2neuron(1, end+1) = v; ui.row2segment(1, end+1) = 0;
        if ~oedb.expandneuron{u}(v) || (oedb.opt.expandonlyselected && (strcmpi(selectiontype,'dataset') || v ~= ui.selectedneuron)), continue; end;
        for w = 1:oedb.nsegments{u}(v)
            listnames{1, end+1} = [tabstr tabstr oedb.segmentnames{u}{v}{w}];
            ui.row2dataset(1, end+1) = u; ui.row2neuron(1, end+1) = v; ui.row2segment(1, end+1) = w;
        end
    end
end
if isempty(listnames)
    listnames = {' '};
    dv = 1;
else
    switch selectiontype
        case 'dataset'
            dv = find(ui.row2dataset == ui.selecteddataset & ui.row2neuron == 0 & ui.row2segment == 0);
        case 'neuron'
            dv = find(ui.row2dataset == ui.selecteddataset & ui.row2neuron == ui.selectedneuron & ui.row2segment == 0);
        case 'segment'
            dv = find(ui.row2dataset == ui.selecteddataset & ui.row2neuron == ui.selectedneuron & ui.row2segment == ui.selectedsegment);
    end
    assert(numel(dv) < 2, 'data indexing error');
    if isempty(dv)
        dv = max(1, min(get(handles.data_listbox,'value'), numel(listnames)));
    end
end
setappdata(handles.oedatabrowser_figure,'ui',ui);
set(handles.data_listbox, 'value', dv, 'string', listnames);

function update_alglist(handles, oedb)
if isempty(oedb.algnames)
    
    v = 1; s = {' '};
    
else
    
    s = oedb.algnames;
    
    if ismember('sbm', s)
        
        v = find(strcmpi('sbm', s));
        assert(numel(v) == 1, 'duplicate algorithm name');
        
    else
        
        v = max(1, min(get(handles.alg_popup, 'value'), numel(s)));
        
    end
    
end
set(handles.alg_popup,'string',s,'value',v);


function ok = update_selection(handles)
ok = false;
[ui, fileinfo, selectiontype, algind, algname, oedb] = oedb_getbasics(handles);
dv = get(handles.data_listbox,'value'); if isempty(dv), return; end;
set([handles.exportdata handles.exportsim],'visible','off');
[ui.selecteddataset, ui.selectedneuron, ui.selectedsegment] = deal(0);
if oedb.algavailable(algind)
    
    set(handles.detect_button, 'enable','on');
    
else
    
    set(handles.detect_button, 'enable','off');
    
end
if oedb.ndatasets > 0
    
    ui.selecteddataset = ui.row2dataset(dv);
    if oedb.nneurons(ui.selecteddataset) > 0
        
        ui.selectedneuron = max(1, ui.row2neuron(dv));
        if oedb.nsegments{ui.selecteddataset}(ui.selectedneuron) > 0
            ui.selectedsegment = max(1, ui.row2segment(dv));
        else
            ui.selectedsegment = 0;
        end
        set(handles.exportdata,'visible','on');
        if oedb.issimulation(ui.selecteddataset)
            set(handles.exportsim,'visible','on');
        end
        
    else
        
        [ui.selectedneuron, ui.selectedsegment] = deal(0);
        
    end
    
end
setappdata(handles.oedatabrowser_figure,'ui',ui);
update_ui_from_algoptsandparams(handles);
ok = true;


function selectiontype = get_selectiontype(handles)
ui = getappdata(handles.oedatabrowser_figure,'ui');
dv = get(handles.data_listbox,'value');
if isempty(ui.row2neuron) || isempty(dv) || ~ui.row2neuron(dv)
    selectiontype = 'dataset';
elseif ~ui.row2segment(dv)
    selectiontype = 'neuron';
else
    selectiontype = 'segment';
end

function varargout = oedatabrowser_OutputFcn(hObject, eventdata, handles)
try
    varargout{1} = handles.output;
    ou = get(handles.oedatabrowser_figure,'units');
    set(handles.oedatabrowser_figure,'units','normalized');
    p = get(handles.oedatabrowser_figure,'position');
    p(1) = max(p(1), 0.01); p(2) = max(p(2), 0.05); p(3) = min(p(3), 0.99 - p(1)); p(4) = min(p(4), 0.9 - p(2));
    set(handles.oedatabrowser_figure,'position',p);
    set(handles.oedatabrowser_figure,'units',ou);
catch %#ok<CTCH>
    varargout{1} = [];
end

function estep_button_Callback(hObject, eventdata, handles)
mi = ui_getval(handles.sbm_iter_edit);
ui_setval(handles.sbm_iter_edit, 0);
detect_button_Callback(handles.detect_button, eventdata, handles);
ui_setval(handles.sbm_iter_edit, mi);


function [ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles)
ui = getappdata(handles.oedatabrowser_figure,'ui');
fileinfo = getappdata(handles.oedatabrowser_figure,'fileinfo');
selectiontype = get_selectiontype(handles);
algind = get(handles.alg_popup,'value');
algname = lower(listboxchoice(handles.alg_popup));
if nargout < 6, return; end
oedb = getappdata(handles.oedatabrowser_figure,'oedb');
segmentlist = [];
if ~(ui.selecteddataset && ui.selectedneuron && ui.selectedsegment), return; end
switch selectiontype
    case 'dataset'
        %segmentlist = [];
    case 'neuron'
        segmentlist = 1:oedb.nsegments{ui.selecteddataset}(ui.selectedneuron);
    case 'segment'
        segmentlist = ui.selectedsegment;
end


function neuronlist = get_selected_neurons(oedb, ui, selectiontype)
if strcmpi(selectiontype, 'dataset')
    
    neuronlist = 1:oedb.nneurons(ui.selecteddataset);
    
else
    
    neuronlist = ui.selectedneuron;
    
end


function [sctrue, scest, tsc, neuronlist, neuron_index] = fetch_selected_trueandest_tandsc(oedb, fileinfo, ui, selectiontype, algind)
neuronlist = get_selected_neurons(oedb, ui, selectiontype);
[sctrue, scest, tsc] = deal({});
neuron_index = [];
for n = neuronlist
    
    segmentlist = get_selected_segments(oedb, ui, selectiontype, n);
    
    oerec = fetch_neurondata(oedb, fileinfo, ui.selecteddataset, n);
    
    for s = segmentlist
        
        data = oerec.data(s);
        results = fetch_results(oedb, fileinfo, ui.selecteddataset, n, s, algind);
        
        %retrieve true and estimated spike counts / time base:
        [sctrue{1, end + 1}, scest{1, end + 1}, tsc{1, end + 1}] = retrieve_spikecounts_for_stats(data, results);
        neuron_index(1, end + 1) = n;
        
    end
    
end


function scalefactorcurves_Callback(hObject, eventdata, handles)
nscalefacs = 100;
minscalefac = 0.05;
maxscalefac = 10;
scalefacs = exp(linspace(log(minscalefac), log(maxscalefac), nscalefacs));
if any(scalefacs > 1) && any(scalefacs < 1) %make it so that one of the scale factors is 1
    
    [~, q] = min(abs(scalefacs - 1));
    scalefacs = scalefacs / scalefacs(q);
    
end

[ui, fileinfo, selectiontype, algind, algname, oedb] = oedb_getbasics(handles);
[sctrue, scest, tsc, neuronlist, neuron_index] = fetch_selected_trueandest_tandsc(oedb, fileinfo, ui, selectiontype, algind);


if ui_ison(handles.scalefactorsdnorm) %scale all neurons' signals so they have the same s.d. value
    
    sdvals = nan(1, numel(neuronlist));
    for n = neuronlist
        
        g = neuron_index == n;
        sdvals(neuronlist == n) = std(cat(1, scest{g}));
        scest(g) = cellfun(@(v) v / sdvals(neuronlist == n), scest(g), 'uniformoutput', false); %divide by this neuron's s.d. value
        
    end
    
    scest = cellfun(@(v) v * mean(sdvals), scest, 'uniformoutput', false); %multiply by the mean s.d. value over all neurons
    
end


[detrate, fprate] = deal(nan(numel(scalefacs), numel(neuronlist)));
wb = waitbar(0, 'computing scaled-data DR/FP curves');
for ii = 1:numel(scalefacs)
    
    scest_scaled = cellfun(@(v) v * scalefacs(ii), scest, 'uniformoutput', false);
    
    [npairs, nspikes_true, nspikes_est, Tseg] = test_spike_recon_fromsc(scest_scaled, sctrue, tsc, oedb.opt.tol_sec);
    
    for n = neuronlist
        
        g = neuron_index == n;
        
        detrate(ii, neuronlist == n) = sum(npairs(g)) / sum(nspikes_true(g));
        
        fprate(ii, neuronlist == n) = sum(nspikes_est(g) - npairs(g)) / sum(Tseg(g));
        
    end
    try 
        waitbar(ii / numel(scalefacs), wb);
    catch
    end
    
end
if isvalid(wb); close(wb); end

f = figure;
ax = axes('parent','f');
plot(fprate, detrate, 'b.-', 'parent', ax);
hold on;
[mfp, mdr] = deal(mean(fprate, 2), mean(detrate, 2));
plot(mfp, mdr, 'k.-', 'linewidth', 2, 'parent', ax);
if any(scalefacs == 1)
    
    hold on;
    plot(mfp(scalefacs == 1), mdr(scalefacs == 1), 'ok', 'parent', ax);
    
end

ylabel('Detection rate');

xlabel('FP [Hz]');


titlestr = sprintf('scalefac curves of FP/DR for %s', algname);

if ui_ison(handles.scalefactorsdnorm)
    
    titlestr = [titlestr ' (normalized)'];
    
end
title(titlestr);
assignin('base','fprate',fprate);
assignin('base','detrate',detrate);


function analyze_roc_Callback(hObject, eventdata, handles)
nrocspts = 100;
[ui, fileinfo, selectiontype, algind, algname, oedb] = oedb_getbasics(handles);

[sctrue, scest, tsc, neuronlist, neuron_index] = fetch_selected_trueandest_tandsc(oedb, fileinfo, ui, selectiontype, algind);

if ui_ison(handles.roc_sdnorm)
    %normalize each neuron's values by the standard deviation
    sdvals = nan(1, numel(neuronlist));
    [minthreshval, maxthreshval] = deal(nan);
    
    for n = neuronlist
        
        g = neuron_index == n;
        scestvals = cat(1, scest{g});
        sdvals(neuronlist == n) = std(scestvals);
        scestvals_norm = scestvals / sdvals(neuronlist == n);
        maxthreshval = max(maxthreshval, max(scestvals_norm));
        minthreshval = min(minthreshval, min(scestvals_norm));
        
    end
    
else
    
    minthreshval = min(cat(1, scest{:}));
    maxthreshval = max(cat(1, scest{:}));
    
end

threshvals = linspace(minthreshval, maxthreshval, nrocspts);

[detrate, fprate] = deal(nan(numel(threshvals), numel(neuronlist)));
wb = waitbar(0, 'computing ROC curve');
for ii = 1:numel(threshvals)
    
    scest_thresholded = cell(1, numel(scest));
    for u = 1:numel(scest)
        
        if ui_ison(handles.roc_sdnorm)
            
            tv = threshvals(ii) * sdvals(neuron_index(u));
            
        else
            
            tv = threshvals(ii);
            
        end
        
        scest_thresholded{u} = scest{u} .* double(scest{u} > tv);
        
    end
    
    if ui_ison(handles.roc_binarize)
        
        scest_thresholded = cellfun(@(v) double(v > 0), scest_thresholded, 'uniformoutput', false);
        sctrue_processed  = cellfun(@(v) double(v > 0),            sctrue, 'uniformoutput', false);
        
    else
        
        sctrue_processed  = sctrue;
        
    end
    [npairs, nspikes_true, nspikes_est, Tseg] = test_spike_recon_fromsc(scest_thresholded, sctrue_processed, tsc, oedb.opt.tol_sec);
    
    for n = neuronlist
        
        g = neuron_index == n;
        
        detrate(ii, neuronlist == n) = sum(npairs(g)) / sum(nspikes_true(g));
        
        if ui_ison(handles.roc_binarize)
            
            ntimepoints = sum(cellfun(@numel, sctrue_processed(g)));
            npositive   = sum(cellfun(@sum,   sctrue_processed(g)));
            fprate(ii, neuronlist == n) = sum(nspikes_est(g) - npairs(g)) / (ntimepoints - npositive);
            
        else
            
            fprate(ii, neuronlist == n) = sum(nspikes_est(g) - npairs(g)) / sum(Tseg(g));
            
        end
        
    end
    try
        waitbar(ii / numel(threshvals), wb);
    catch
    end
    
end
if isvalid(wb); close(wb); end

figure;
plot(fprate, detrate, 'b.-');
hold on;
plot(mean(fprate, 2), mean(detrate, 2), 'k.-', 'linewidth', 2);

ylabel('Detection rate');
if ui_ison(handles.roc_binarize)
    
    xlabel('FP rate');
    
else
    
    xlabel('FP [Hz]');
    
end
titlestr = sprintf('ROC for %s', algname);

if ui_ison(handles.roc_binarize) && ui_ison(handles.roc_sdnorm)
    
    titlestr = [titlestr ' (normalized, binarized)'];
    
elseif ui_ison(handles.roc_binarize)
    
    titlestr = [titlestr ' (binarized)'];
    
elseif ui_ison(handles.roc_sdnorm)
    
    titlestr = [titlestr ' (normalized)'];
    
end
title(titlestr);
assignin('base','fprate',fprate);
assignin('base','detrate',detrate);


function detect_button_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, algind, algname, oedb] = oedb_getbasics(handles);
if ~(ui.selecteddataset && ui.selectedneuron && ui.selectedsegment), return; end

[saveafter, saveafterevery] = deal(false);
if strcmpi(selectiontype, 'dataset')
    
    answer = questdlg('Save when finished?', 'AP detection', 'At end', 'After each neuron', 'No', 'No');
    if isempty(answer), return; end
    
    if ~strcmpi(answer, 'No')
        
        saveafter = true;
        saveafterevery = strcmpi(answer, 'After each neuron');
        [f, p] = putfile_fromp(fileinfo.odbdir, '*.odb', 'Save optical/electrical database', fileinfo.odbfile);
        if isnumeric(f), return; end
        
    end
    
    prev_results_present = false(1, oedb.nneurons(ui.selecteddataset));
    for neuronind = 1:oedb.nneurons(ui.selecteddataset)
        
        for s = 1:oedb.nsegments{ui.selecteddataset}(neuronind)
            
            r = fetch_results(oedb, fileinfo, ui.selecteddataset, neuronind, s, algind);
            if ~isempty(r.spikecounts) || ~any(isnan(r.spiketimes_window))
                
                prev_results_present(neuronind) = true;
                break;
                
            end
            
        end
        
    end
    neuronlist = 1:oedb.nneurons(ui.selecteddataset);
    if any(prev_results_present)
        
        answer = questdlg(sprintf('Overwrite previous results in %d neurons?', sum(prev_results_present)), 'AP detection', 'Overwrite', 'Skip', 'Overwrite');
        if isempty(answer), return; end
        if strcmpi(answer, 'Skip')
            
            neuronlist(prev_results_present) = [];
            
        end
        if isempty(neuronlist), return; end
        
    end
    
end

oedb.modified = true;
set(handles.oedatabrowser_figure,'pointer','watch'); drawnow;
if ui_ison(handles.opt_profile), profile on; end

switch selectiontype
    
    case 'dataset'
        
        wb = waitbar(0, ['Analyzing neuron 1 / ' num2str(numel(neuronlist))]); drawnow;
        for j= 1:numel(neuronlist)
            try
                waitbar(j / numel(neuronlist), wb, ['Analyzing neuron ' num2str(j) ' / ' num2str(numel(neuronlist))]); drawnow;
            catch
                wb = waitbar(j / numel(neuronlist), ['Analyzing neuron 1 / ' num2str(numel(neuronlist))]); drawnow;
            end
            neuronind = neuronlist(j);
            
            oedb = oedatabrowser_detectspikes(oedb, handles, fileinfo, algind, algname, ui.selecteddataset, neuronind); %all segments
            pause(0.05);
            
            if saveafterevery
                
                setappdata(handles.oedatabrowser_figure,'oedb',oedb);
                [oedb, fileinfo] = save_fromfig(handles.oedatabrowser_figure, p, f, oedb, fileinfo);
                
            end
            
        end
        if isvalid(wb); close(wb); end
        
    case 'neuron'
        
        oedb = oedatabrowser_detectspikes(oedb, handles, fileinfo, algind, algname, ui.selecteddataset, ui.selectedneuron); %all segments
        
    case 'segment'
        
        oedb = oedatabrowser_detectspikes(oedb, handles, fileinfo, algind, algname, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment);
        
end

setappdata(handles.oedatabrowser_figure,'oedb',oedb);
if saveafter && ~saveafterevery
    
    [oedb, fileinfo] = save_fromfig(handles.oedatabrowser_figure, p, f, oedb, fileinfo);
    
end
if ui_ison(handles.opt_profile), profile off; profile report; end
update_ui_from_algoptsandparams(handles);
show_data(handles);
set(handles.oedatabrowser_figure,'pointer','arrow');


function clear_button_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, algind, ~, oedb] = oedb_getbasics(handles);
if ~(ui.selecteddataset && ui.selectedneuron && ui.selectedsegment), return; end
switch selectiontype
    case 'dataset'
        oedb = oedb_clear_results(oedb, fileinfo, handles, algind, ui.selecteddataset);
    case 'neuron'
        oedb = oedb_clear_results(oedb, fileinfo, handles, algind, ui.selecteddataset, ui.selectedneuron);
    case 'segment'
        oedb = oedb_clear_results(oedb, fileinfo, handles, algind, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment);
end
oedb.modified = true;
setappdata(handles.oedatabrowser_figure,'oedb',oedb);
update_ui_from_algoptsandparams(handles);
show_data(handles);
update_figurename(handles.oedatabrowser_figure, oedb, fileinfo);


function fn = Pfields_ui
fn = {'gain' 'zeta' 'sigma_r' 'fr' 'R' 'S'};


function clear_Pui(handles)
Pfields_toui = Pfields_ui;
for fi = 1:numel(Pfields_toui)
    fn = Pfields_toui{fi};
    ui_setval([handles.([fn '_edit']) handles.([fn '_text'])], NaN);
end


function vok = sanitize_param_value(algname, param_name, v)
vok = NaN;
if ~strcmpi(algname, 'sbm'), return; end; %FIXME

v = real(v);
if isinf(v), v = NaN; end
switch param_name
    
    case {'sigma_r' 'fr' 'R'}
        
        v = max(v, 1e-4); %should probably have some more explicit minima due to numerical stability etc. FIXME
        
    case 'gain'
        
        v = max(v, 0);
        
    case 'S'
        
        v = min(max(v, 0.01), 100);
        
    case 'zeta'
        
        v = max(v, eps);
        
    otherwise
        
        warning('oedatabrowser:unrecognizedparameter','unrecognized param'); return;
        
end
vok = v;


function Pfromui(handles, hObject)
[ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles);

if ~strcmpi(algname, 'sbm')
    
    set(hObject, 'string', '');
    return; %FIXME
    
end
fn = get(hObject,'tag');
assert(strcmpi(fn(end-4:end),'_edit'), 'invalid ui object');
fn = fn(1:end-5);

Pv = ui_getval(hObject);
assert(ismember(fn, Pfields_ui));
Pv = sanitize_param_value(algname, fn, Pv);
if isnan(Pv)
    
    update_ui_from_algoptsandparams(handles);
    return;
    
end

if strcmpi(selectiontype, 'dataset')
    
    neuronlist = 1:oedb.nneurons(ui.selecteddataset);
    
else
    
    neuronlist = ui.selectedneuron;
    
end

for n = neuronlist
    
    if strcmpi(selectiontype, 'dataset')
        
        segmentlist = 1:oedb.nsegments{ui.selecteddataset}(n);
        
    end
    
    for s = segmentlist
        
        r = fetch_results(oedb, fileinfo, ui.selecteddataset, n, s, algind);
        r.params.(fn) = Pv;
        oedb = oedb_assign_results(oedb, fileinfo, ui.selecteddataset, n, s, algind, r);
        
    end
    
end

oedb.modified = true;
setappdata(handles.oedatabrowser_figure,'oedb',oedb);
show_data(handles);


function h = simonly_uiobj(handles)
h = [handles.gain_text handles.S_text handles.R_text handles.zeta_text ...
    handles.sigma_r_text handles.fr_text ...
    ];


function s = fieldvalstring_fromcellofstructs(fn, c)
s = ''; %default, no values

fieldok = cellfun(@(v) isa(v, 'struct') && isfield(v, fn) && (isa(v.(fn), 'string') || ~any(isnan(v.(fn)))), c);
if ~any(fieldok), return; end
c = c(fieldok);

cisscalar = cellfun(@(v) isscalar(v.(fn)), c);
if ~all(cisscalar)
    
    s = '<array>';
    return;
    
end

x = cellfun(@(v) double(v.(fn)), c); %field values
if all(x == x(1)) && all(fieldok) %if some values are missing, show as 'mixed'
    
    s = sprintf('%7.5g', x(1));
    
else
    
    s = 'mixed';
    
end


function update_ui_from_algoptsandparams(handles)
[ui, fileinfo, selectiontype, algind, algname, oedb] = oedb_getbasics(handles);

ui_setval(handles.LL_logmarglik_text, NaN);
clear_Pui(handles);

Pfields_toui = Pfields_ui;
h_simonly = simonly_uiobj(handles);
h_sbmonly = [handles.estep_button handles.parameterestimation_popup handles.parameterestimation_text handles.nparticles_edit ...
    handles.nparticles_text handles.sbm_iter_edit handles.sbm_iter_text handles.maxjitter_edit handles.maxjitter_text];

set([h_simonly h_sbmonly],'visible','off');

if ~all([ui.selecteddataset ui.selectedneuron ui.selectedsegment]), return; end

oerec = fetch_neurondata(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron);

ui_setval(handles.settingsindex_edit,   oerec.data(ui.selectedsegment).info.settingsindex);
ui_setval(handles.focalplaneindex_edit, oerec.data(ui.selectedsegment).info.focalplaneindex);

neuronlist = get_selected_neurons(oedb, ui, selectiontype);

[params, opts, params_sim] = deal({});
log_marglik = nan(0, 1);

for n = neuronlist
    
    segmentlist = get_selected_segments(oedb, ui, selectiontype, n);
    
    %get results of algorithm run etc.
    for s = 1:numel(segmentlist)
        
        results(s) = fetch_results(oedb, fileinfo, ui.selecteddataset, n, segmentlist(s), algind);
        
    end
    
    %aggregate in cell arrays instead of struct arrays since some fields may be missing for some data
    params = cat(2, params, {results.params});
    opts   = cat(2,   opts, {results.opts});
    
    %get cost function
    for s = 1:numel(segmentlist)
        
        if isfield(results(s).outputvars, 'lik')
            
            LL = results(s).outputvars.lik.logmarglik;
            
        else
            
            LL = NaN;
            
        end
        
        log_marglik(1, end + 1) = LL;
        
    end
    
    %get simulation info
    if oedb.issimulation(ui.selecteddataset)
        
        siminfo = fetch_siminfo(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron);
        %aggregate in cell array instead of struct array since some fields may be missing for some data
        params_sim = cat(2, params_sim, {siminfo(segmentlist).P});
        
    end
    
end

for fi = 1:numel(Pfields_toui)
    
    fn = Pfields_toui{fi};
    
    set(handles.([fn '_edit']), 'string', fieldvalstring_fromcellofstructs(fn, params));
    set(handles.([fn '_text']), 'string', fieldvalstring_fromcellofstructs(fn, params_sim));
    
end

if oedb.issimulation(ui.selecteddataset)
    
    set(h_simonly, 'visible', 'on');
    
end

if strcmpi(algname, 'sbm')
    
    set(h_sbmonly, 'visible', 'on');
    
    optsui = opts{1};
    
    if isfield(optsui, 'shotnoise') %otherwise, leave the gui in whatever state it was
        
        if optsui.shotnoise
            
            set(handles.sbm_shotnoise,'checked','on');
            
        else
            
            set(handles.sbm_shotnoise,'checked','off');
            
        end
        
    end
    
end

if all(isnan(LL))
    
    ui_setval(handles.LL_logmarglik_text, NaN);
    
else
    
    ui_setval(handles.LL_logmarglik_text, sum(LL(~isnan(LL))));
    
end


function opts = sbm_optsfromui(handles)
opts = orderfields(struct( ...
    'min_iter_paramest', ui_getval(handles.sbm_iter_edit), 'Nparticles', ui_getval(handles.nparticles_edit) ...
    ,'filtersmoother_window', ui_getval(handles.fswin_edit) ...
    ,'verbose', verbositylevel_query(handles) ...
    ,'shotnoise', ui_ison(handles.sbm_shotnoise) ...
    ,'Nspikehist',0 ...
    ,'parameterestimation', lower(listboxchoice(handles.parameterestimation_popup)) ...
    ,'usegpu',ui_ison(handles.sbm_usegpu) ...
    ,'gfitst',ui_ison(handles.sbm_moments2jitter) ...
    ,'profile_firing_rate', ui_ison(handles.sbm_profileFR) ...
    ,'profile_zeta', ui_ison(handles.sbm_profilezeta) ...
    ,'maxjitter_ms', ui_getval(handles.maxjitter_edit) ...
    ,'resample_to_proposal', true ...
    ,'multiroundpf',true ...
    ));


function oedb = oedatabrowser_detectspikes(oedb, handles, fileinfo, algind, algname, datasetind, neuronind, segmentlist)
if ~exist('segmentlist','var')
    
    segmentlist = 1:oedb.nsegments{datasetind}(neuronind);
    
end
if isempty(segmentlist) || ~datasetind || ~neuronind, return; end

oerec = fetch_neurondata(oedb, fileinfo, datasetind, neuronind);
oerec.data = oerec.data(segmentlist);

if ui_ison(handles.onscreenonly) && numel(segmentlist) == 1 %only process onscreen data
    
    [~,it,dt] = extract_data_fromoerec(oerec, 1);
    XL = get(handles.f_axes,'xlim');
    zoommask = reshape(it{1} > XL(1) - dt & it{1} < XL(2) + dt,1,[]);
    oerec.data = oerecdata_extract_framerange(oerec.data, zoommask);
    
end

for s = 1:numel(segmentlist)
    
    prevresults(s) = fetch_results(oedb, fileinfo, datasetind, neuronind, segmentlist(s), algind);
    
end
params = {prevresults.params};
opts = optsfromui(handles);

%------------>> run the algorithm <<------------
results = oedb.algfunctions{algind}(oerec, {params}, opts);
results = results{1};

for segmentind = 1:numel(segmentlist)
    
    oedb = oedb_assign_results(oedb, fileinfo, datasetind, neuronind, segmentlist(segmentind), algind, results(segmentind));
    
end

%update stats:
oedb = oedb_calc_neuron_stats(oedb, fileinfo, algind, datasetind, neuronind);


function opts = optsfromui(handles)
[ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles);
if strcmpi(algname, 'sbm')
    
    opts = sbm_optsfromui(handles);
    
else
    
    opts = struct();
    
end


function exportdata_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, ~, ~, oedb] = oedb_getbasics(handles);
if ~(ui.selecteddataset && ui.selectedneuron && ui.selectedsegment), return; end
switch selectiontype
    case 'dataset'
        oerec = fetch_dataset_oerecarray(oedb, fileinfo, ui.selecteddataset); %#ok<NASGU>
        defaultname = ['dataset_' oedb.datasetnames{ui.selecteddataset} '.mat'];
    case 'neuron'
        oerec = fetch_neurondata(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron); %#ok<NASGU>
        defaultname = ['dataset_' oedb.datasetnames{ui.selecteddataset} '_neuron_' num2str(ui.selectedneuron) '.mat'];
    case 'segment'
        oerec = fetch_neurondata(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron);
        oerec.data = oerec.data(ui.selectedsegment);
        defaultname = ['dataset_' oedb.datasetnames{ui.selecteddataset} '_neuron_' num2str(ui.selectedneuron) '_segment_' num2str(ui.selectedsegment) '.mat'];
end
[f,p] = putfile_fromp(fileinfo.odbdir, '*.mat', ['Save optical/electrical data (oerec) for ' selectiontype], defaultname);
if isnumeric(f), return; end;
fileinfo.odbdir = p;
set(handles.oedatabrowser_figure,'pointer','watch'); drawnow;
save([p f],'oerec','-mat');
set(handles.oedatabrowser_figure,'pointer','arrow');


function file_save_Callback(hObject, eventdata, handles)
fileinfo = getappdata(handles.oedatabrowser_figure,'fileinfo');
[f ,p] = putfile_fromp(fileinfo.odbdir, '*.odb', 'Save optical/electrical database', fileinfo.odbfile);
if isnumeric(f), return; end
set(handles.oedatabrowser_figure,'pointer','watch'); drawnow;
save_fromfig(handles.oedatabrowser_figure, p, f);
set(handles.oedatabrowser_figure, 'pointer','arrow'); drawnow;


function [oedb, fileinfo] = save_fromfig(figh, p, f, oedb, fileinfo)
if nargin < 4
    oedb = getappdata(figh, 'oedb');
end
if nargin < 5
    fileinfo = getappdata(figh, 'fileinfo');
end
ui = getappdata(figh, 'ui');
[oedb, fileinfo] = oedb_save(p, f, oedb, fileinfo, ui);
setappdata(figh, 'oedb', oedb);
setappdata(figh, 'fileinfo', fileinfo);
update_figurename(figh, oedb, fileinfo);


function file_load_Callback(hObject, eventdata, handles)
oedb = getappdata(handles.oedatabrowser_figure, 'oedb');
fileinfo = getappdata(handles.oedatabrowser_figure, 'fileinfo');
[f, p] = getfile_fromp(fileinfo.odbdir, '*.odb', 'Load optical/electrical database');
if isnumeric(f), return; end;
set(handles.oedatabrowser_figure,'pointer','watch'); drawnow;

prevalgs = orderfields(struct('algfunctions',{oedb.algfunctions},'alginfo',oedb.alginfo,'algnames',{oedb.algnames},'trainfunctions',{oedb.trainfunctions}));

loadok = false;
try
    
    [oedb, fileinfo, ui] = oedb_load([p f], fileinfo, oedb.opt.figsavedir);
    loadok = true;
    
catch ex
    
    uiwait(errordlg(ex.message, 'Failed to load file', 'modal'));
    
end
set(handles.oedatabrowser_figure,'pointer','arrow');
if ~loadok, return; end

ui.row2neuron = []; %otherwise we can get errors. this will be redone by update_ui

%see which algorithms from the file are available in the present context
nalg = numel(oedb.alginfo);
oedb.algavailable = false(1, nalg);
for algind = 1:nalg
    
    matchind = find(strcmpi(prevalgs.algnames, oedb.algnames{algind}),1);
    if isempty(matchind), continue; end;
    [oedb.algfunctions{algind}, oedb.alginfo(algind), oedb.trainfunctions{algind}] = deal(prevalgs.algfunctions{matchind}, prevalgs.alginfo(matchind), prevalgs.trainfunctions{matchind});
    oedb.algavailable(algind) = true;
    
end

%add any algorithms available on this system that were not mentioned in the file

for prevalgind = 1:numel(prevalgs.alginfo)
    
    if any(strcmpi(oedb.algnames, prevalgs.algnames{prevalgind})), continue; end; %algorithm was already present. FIXME check versions?
    
    oedb.alginfo(1, end+1) = prevalgs.alginfo(prevalgind);
    oedb.algnames{1, end+1} = prevalgs.algnames{prevalgind};
    oedb.algfunctions{1, end + 1} = prevalgs.algfunctions{prevalgind};
    oedb.trainfunctions{1, end + 1} = prevalgs.trainfunctions{prevalgind};
    oedb.algavailable(1, end+1) = true;
    
    oedb.results(end + 1, 1:oedb.ndatasets) = cell(1, oedb.ndatasets);
    oedb.stats.bydataset(end + 1, :) = repmat(empty_statstruct, 1, oedb.ndatasets);
    
    for u = 1:oedb.ndatasets
        
        if oedb.opt.storedatainmemory
            
            oedb.results{end, u} = apdet_resultsstruct(oedb.nsegments{u});
            
        else
            
            oedb.results{end, u} = cell(1, oedb.nneurons(u));
            
        end
        
        oedb.stats.byneuron{u}(end + 1, :) = repmat(empty_statstruct, 1, oedb.nneurons(u));
        
        for n = 1:oedb.nneurons(u)
            
            oedb.stats.bysegment{u}{n}(end + 1, :) = repmat(empty_statstruct, 1, oedb.nsegments{u}(n));
            
        end
        
    end
    
end

%store results:
setappdata(handles.oedatabrowser_figure,'oedb',oedb);
setappdata(handles.oedatabrowser_figure,'fileinfo',fileinfo);
setappdata(handles.oedatabrowser_figure,'ui',ui);
%update gui:
update_ui(handles, oedb, fileinfo);
show_data(handles);
zoomoutfull(handles);


function oedatabrowser_figure_WindowScrollWheelFcn(hObject, eventdata, handles)
oedb_swf(hObject, eventdata, handles);


function oedb_swf(hObject, eventdata, handles)
%scroll-wheel function
zoomax = [handles.f_axes handles.spike_axes handles.ca_axes];
ui = getappdata(handles.oedatabrowser_figure,'ui');
if any([ui.selecteddataset ui.selectedneuron ui.selectedsegment] == 0), return; end;
[clickax, clickdatapos, ~, ~, ~] = ui_axes_clickinfo(hObject);
if ismember(clickax, zoomax)
    
    oedb = getappdata(hObject,'oedb');
    fileinfo = getappdata(handles.oedatabrowser_figure,'fileinfo');
    [~,it,~,~] = oedb_fetch_timeseries(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment); %f,it column vectors; ns row vector
    dt = median(diff(it));
    prevXL = get(clickax,'xlim');
    if eventdata.VerticalScrollCount < 0 %zoom in
        newXL = clickdatapos(1) + (prevXL - clickdatapos(1)) / 2;
        if diff(newXL) < 25 * dt, return; end
    else %zoom out
        newXL = clickdatapos(1) + (prevXL - clickdatapos(1)) * 2;
    end
    if newXL(1) < it(1) - 2 * dt
        newXL = newXL + (it(1) - 2 * dt - newXL(1));
    end
    if newXL(2) > it(end) + 2 * dt
        newXL = newXL - (newXL(2) - it(end) - 2 * dt);
    end
    newXL = max(it(1) - 2 * dt, min(it(end) + 2 * dt, newXL));
    set(zoomax,'xlim',newXL);
    
end
drawnow expose;
pause(0.025);


function oedatabrowser_figure_WindowButtonDownFcn(hObject, eventdata, handles)
ui = getappdata(handles.oedatabrowser_figure,'ui');
if any([ui.selecteddataset ui.selectedneuron ui.selectedsegment] == 0), return; end;
zoomax = [handles.f_axes handles.spike_axes handles.ca_axes];
[clickax, clickdatapos, Lclick, Rclick, Oclick] = ui_axes_clickinfo(hObject);
if ismember(clickax, zoomax)
    
    oedb = getappdata(hObject,'oedb');
    fileinfo = getappdata(handles.oedatabrowser_figure,'fileinfo');
    [~,it,~,~] = oedb_fetch_timeseries(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment); %f,it column vectors; ns row vector
    dt = median(diff(it));
    if Rclick || true
        rbbox;
        cp = get(clickax,'currentpoint');
        newXL = sort([clickdatapos(1) cp(1)]);
        prevXL = get(clickax,'xlim');
        if diff(newXL) / diff(prevXL) < 0.02
            newXL = mean(prevXL) + (prevXL - mean(prevXL)) * 2; %zoom out
        end
        if diff(newXL) < 10 * dt
            newXL = mean(newXL) + [-1 1] * 5 * dt;
        end
        newXL = min(max(newXL, it(1) - 2 * dt), it(end) + 2 * dt);
        if diff(newXL) <= 0
            
            newXL = [it(1) - 2 * dt, it(end) + 2 * dt];
            
        end
        set(zoomax,'xlim',newXL);
        
    elseif Oclick
        
        zoomoutfull(handles);
        
    end
    
end


function oedatabrowser_figure_ResizeFcn(hObject, eventdata, handles)


function data_listbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function data_listbox_Callback(hObject, eventdata, handles)
if ~update_selection(handles), return; end; %not sure why this callback gets called sometimes while loading a file
oedb = getappdata(handles.oedatabrowser_figure,'oedb');
if oedb.opt.expandonlyselected
    update_datalist(handles, oedb);
end
show_data(handles);
zoomoutfull(handles);


function alg_popup_Callback(hObject, eventdata, handles)
update_selection(handles);
show_data(handles);


function alg_popup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function cheat_checkbox_Callback(hObject, eventdata, handles)
show_data(handles);


function sbm_iter_edit_Callback(hObject, eventdata, handles)


function sbm_iter_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nparticles_edit_Callback(hObject, eventdata, handles)

function nparticles_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function oedatabrowser_figure_CloseRequestFcn(hObject, eventdata, handles)
fileinfo = getappdata(handles.oedatabrowser_figure,'fileinfo');
try
    oedb_cleartmpdir(fileinfo);
    rmdir(fileinfo.tmpdir);
catch ex
    uiwait(errordlg(ex.message, 'Cannot clear/remove tmp directory', 'modal'));
end
delete(hObject);


function opt_storedatainmemory_Callback(hObject, eventdata, handles)
set(handles.oedatabrowser_figure,'pointer','watch'); drawnow;
oedb = getappdata(handles.oedatabrowser_figure,'oedb');
fileinfo = getappdata(handles.oedatabrowser_figure,'fileinfo');
try
    if oedb.opt.storedatainmemory
        
        oedb = oedb_movedatatodisk(oedb, fileinfo);
        
    else
        
        oedb = oedb_movedatatomemory(oedb, fileinfo);
        
    end
catch ex
    
    uiwait(errordlg(ex.message, 'Failed to move data', 'modal'));
    set(handles.oedatabrowser_figure,'pointer','arrow');
    return;
    
end
oedb.modified = true;
setappdata(handles.oedatabrowser_figure, 'oedb',oedb);
update_ui_oedbopts(handles, oedb);
set(handles.oedatabrowser_figure,'pointer','arrow');


function fr_edit_Callback(hObject, eventdata, handles) %#ok<*INUSL>
[ui, fileinfo, ~, ~, algname, oedb] = oedb_getbasics(handles);
if ~strcmpi(algname, 'sbm'), return; end
if ~all([ui.selecteddataset ui.selectedneuron ui.selectedsegment]), ui_setval(hObject,NaN); return; end;
%ensure that spiking rate parameter is not less than one spike per entire recording length for this neuron
npts = 0;
for segment = 1:oedb.nsegments{ui.selecteddataset}(ui.selectedneuron)
    [~,it,~,~] = oedb_fetch_timeseries(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, segment); %f,it column vectors; ns row vector
    npts = npts + numel(it);
end
if npts == 0, ui_setval(hObject,NaN); return; end;
ui_setval(hObject, max(ui_getval(hObject), 1 / npts));
Pfromui(handles, hObject);


function R_edit_Callback(hObject, eventdata, handles)
Pfromui(handles, hObject);

function S_edit_Callback(hObject, eventdata, handles)
Pfromui(handles, hObject);

function gain_edit_Callback(hObject, eventdata, handles)
Pfromui(handles, hObject);

function gain_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function zeta_edit_Callback(hObject, eventdata, handles)
Pfromui(handles, hObject);

function zeta_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function sigma_r_edit_Callback(hObject, eventdata, handles)
Pfromui(handles, hObject);


function sigma_r_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function opt_profile_Callback(hObject, eventdata, handles)
if ui_ison(hObject)
    set(hObject,'checked','off');
else
    set(hObject,'checked','on');
end


function fr_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function exponlyselected_Callback(hObject, eventdata, handles)
oedb = getappdata(handles.oedatabrowser_figure,'oedb');
if ui_ison(hObject)
    set(hObject,'checked','off');
else
    set(hObject,'checked','on');
end
oedb.opt.expandonlyselected = ui_ison(hObject);
setappdata(handles.oedatabrowser_figure,'oedb',oedb);
update_datalist(handles, oedb);


function file_makevars_Callback(hObject, eventdata, handles)
dataout(handles);


function oerec = neurondatafunc_external(figh, dataset, neuron)
handles = guidata(figh);
[ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles);
oerec = fetch_neurondata(oedb, fileinfo, dataset, neuron);

function dataout(handles)
[ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles);
results = fetch_results(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment, algind);
[f,it,dt,ns,indicatorstring,nA2D,st,f_removed] = oedb_fetch_timeseries(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment); %f,it column vectors; ns row vector

assignin('base','figh',handles.oedatabrowser_figure); assignin('base','ui',ui); assignin('base','algind',algind);
assignin('base','dt',dt); assignin('base','results',results); assignin('base', 'f', f); assignin('base', 'f_removed', f_removed); assignin('base','ns',ns); assignin('base','it',it); assignin('base','st',st); assignin('base','nA2D',nA2D);
assignin('base', 'neurondatafunc', @neurondatafunc_external);

function tolvsaccuracy_Callback(hObject, eventdata, handles)
answer = inputdlg({'min. tol' 'max. tol' 'number of steps'}, 'Correlation vs. timing tolerance', 1, {'0.02' '0.3' '29'}); if isempty(answer), return; end
tol_list = linspace(str2double(answer{1}), str2double(answer{2}), str2double(answer{3}));
[csc, cst, drsc, fpsc, drst, fpst] = deal(nan(size(tol_list)));

[ui, fileinfo, selectiontype, algind, algname, oedb] = oedb_getbasics(handles);
neuronlist = get_selected_neurons(oedb, ui, selectiontype);
ds = ui.selecteddataset;

%fixme: option to do multiple datasets / algorithms at once
wb = waitbar(0, 'Calculating correlation, detection rate and false postive rate vs. timing tolerance');
for k = 1:numel(tol_list)
    
    oedb.opt.tol_sec = tol_list(k);
    
    for j = 1:numel(neuronlist)
        
        n = neuronlist(j);
        
        oedb = oedb_calc_neuron_stats(oedb, fileinfo, algind, ds, n); %fixme: for a single segment, only calculate stats for that segment.
        
    end
    
    if strcmpi(selectiontype, 'segment')
        
        s = oedb.stats.bysegment{ds}{ui.selectedneuron}(algind, ui.selectedsegment);
        
    elseif strcmpi(selectiontype, 'neuron')
        
        s = oedb.stats.byneuron{ds}(algind, ui.selectedneuron);
        
    else %dataset
        
        s = oedb.stats.bydataset(algind, ds);
        
    end
    
    csc(k) = s.corr_spikecounts;
    cst(k) = s.corr_spiketimes;
    drsc(k) = s.det_spikecounts;
    drst(k) = s.det_spiketimes;
    fpsc(k) = s.fp_spikecounts;
    fpst(k) = s.fp_spiketimes;
    
    try
        waitbar(k / numel(tol_list), wb);
    catch
    end
    
end
if isvalid(wb); close(wb); end

figure('name',sprintf('Correlation vs. timing tolerance for %s', algname));
plot(tol_list * 1000, csc, 'k.-'); hold on;
plot(tol_list * 1000, cst, 'm.-');
ylabel('Correlation');
xlabel('Timing tolerance [ms]');

figure('name',sprintf('DR/FP vs. timing tolerance for %s', algname));
subplot(2,1,1);
plot(tol_list * 1000, drsc * 100, 'k.-'); hold on;
plot(tol_list * 1000, drst * 100, 'm.-');
ylabel('Detection rate');
subplot(2,1,2);
plot(tol_list * 1000, fpsc, 'k.-'); hold on;
plot(tol_list * 1000, fpst, 'm.-');
ylabel('False positive rate [Hz]');
xlabel('Timing tolerance [ms]');

assignin('base', 'tol_list', tol_list);
assignin('base', 'csc', csc);
assignin('base', 'cst', cst);
%don't call setappdata on oedb -- these results are to be discarded!


function tol_edit_Callback(hObject, eventdata, handles)
oedb = getappdata(handles.oedatabrowser_figure,'oedb');
fileinfo = getappdata(handles.oedatabrowser_figure,'fileinfo');
tol_sec = ui_getval(hObject) / 1000;
if ~isfinite(tol_sec) || tol_sec <= 0 || imag(tol_sec) ~= 0
    
    ui_setval(hObject, oedb.opt.tol_sec * 1000);
    return;
    
end
oedb.opt.tol_sec = tol_sec;

%redo all stats
for algind = 1:numel(oedb.alginfo)
    for dataset = 1:oedb.ndatasets
        for neuron = 1:oedb.nneurons(dataset)
            
            oedb = oedb_calc_neuron_stats(oedb, fileinfo, algind, dataset, neuron);
            
        end
    end
end
setappdata(handles.oedatabrowser_figure, 'oedb', oedb);
show_data(handles);


function tol_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function sbm_verbosity_Callback(hObject, eventdata, handles)


function verbositylevel_Callback(hObject, eventdata)
set(get(get(hObject,'parent'),'children'),'checked','off');
set(hObject,'checked','on');


function v = verbositylevel_query(handles)
vmenu = get(handles.sbm_verbosity,'children');
vmenu = vmenu(find(ui_ison(vmenu),1));
assert(numel(vmenu) == 1,'GUI state error');
L = get(vmenu,'label');
firstnondigit = find(~ismember(L, num2str((0:9)')'),1);
if isempty(firstnondigit)
    vstr = L;
else
    vstr = L(1:firstnondigit - 1);
end
assert(~isempty(vstr),'GUI state error');
v = str2double(vstr);


function fswin_edit_Callback(hObject, eventdata, handles)


function fswin_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function oedatabrowser_figure_WindowButtonMotionFcn(hObject, eventdata, handles)


function H_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function nonlinfit_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, algind, ~, oedb] = oedb_getbasics(handles);
if any(~[ui.selecteddataset, ui.selectedneuron, ui.selectedsegment]), return; end %no data
if ~strcmpi(oedb.algnames{algind}, 'sbm')
    
    uiwait(errordlg('SBM only', 'Error', 'modal')); return;
    
end

saveafter = false;
if strcmpi(selectiontype, 'dataset')
    
    answer = questdlg('Save odb file when finished?', 'Model fit', 'Yes','No','Cancel','Yes');
    if strcmp(answer, 'Cancel'), return; end
    
    if strcmp(answer, 'Yes')
        
        saveafter = true;
        [f, p] = putfile_fromp(fileinfo.odbdir, '*.odb', 'Save optical/electrical database', fileinfo.odbfile);
        if isnumeric(f), return; end
        
    end
    
end

set(handles.oedatabrowser_figure,'pointer','watch'); drawnow;
if ui_ison(handles.opt_profile), profile on; end

[oedb, canceled] = oedb_nonlinfit(oedb, ui, fileinfo, selectiontype, algind, handles);

if ui_ison(handles.opt_profile), profile off; profile report; end

setappdata(handles.oedatabrowser_figure,'oedb',oedb);
update_ui_from_algoptsandparams(handles);
if saveafter && ~canceled
    
    save_fromfig(handles.oedatabrowser_figure, p, f);
    
end
set(handles.oedatabrowser_figure,'pointer','arrow'); drawnow;


function calc_averages_Callback(hObject, eventdata, handles)
if ui_ison(handles.opt_profile), profile on; end;
oedb = getappdata(handles.oedatabrowser_figure, 'oedb');
fileinfo = getappdata(handles.oedatabrowser_figure,'fileinfo');
oedb_doavcalc(oedb, fileinfo);
if ui_ison(handles.opt_profile), profile off; profile report; end;


function oedatabrowser_figure_WindowKeyPressFcn(hObject, eventdata, handles)
switch eventdata.Key
    case 'd'
        if ismember('control', eventdata.Modifier)
            dataout(handles);
        end
        %     case 'downarrow'
        %         if gco == handles.data_listbox, return; end %otherwise we'd change it twice
        %         set(handles.data_listbox,'value',max(1, get(handles.data_listbox,'value') + 1));
        %         data_listbox_Callback(handles.data_listbox, [], handles);
        %     case 'uparrow'
        %         if gco == handles.data_listbox, return; end %otherwise we'd change it twice
        %         dlstring = get(handles.data_listbox,'string'); if ~isa(dlstring,'cell'), return; end
        %         set(handles.data_listbox,'value',min(numel(dlstring), get(handles.data_listbox,'value') - 1));
        %         data_listbox_Callback(handles.data_listbox, [], handles);
end


function parameterestimation_popup_Callback(hObject, eventdata, handles)

function parameterestimation_popup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function useparallel_Callback(hObject, eventdata, handles)
ui_toggleon(hObject);

function sbm_usegpu_Callback(hObject, eventdata, handles)
ui_toggleon(hObject);

function sbm_moments2jitter_Callback(hObject, eventdata, handles)
ui_toggleon(hObject);

function S_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function R_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function onscreenonly_Callback(hObject, eventdata, handles)
ui_toggleon(hObject);


function sbm_resetgpu_Callback(hObject, eventdata, handles)
cuda5s_devicereset;


function setdatainfofield(fieldname, edith, handles)
[ui, fileinfo, selectiontype, algind, algname, oedb] = oedb_getbasics(handles);
if ~strcmpi(algname, 'sbm') || ~all([ui.selecteddataset ui.selectedneuron ui.selectedsegment]), return; end;
oerec = fetch_dataset_oerecarray(oedb, fileinfo, ui.selecteddataset);
v = ui_getval(edith);
switch selectiontype
    case 'dataset'
        
        for neuronind = 1:oedb.nneurons(ui.selecteddataset)
            
            for segmentind = 1:oedb.nsegments{ui.selecteddataset}(neuronind)
                
                oerec(neuronind).data(segmentind).info.(fieldname) = v;
                
            end
            oedb = modify_neurondata(oedb, fileinfo, ui.selecteddataset, neuronind, oerec(neuronind));
            
        end
        
    case 'neuron'
        
        for segmentind = 1:oedb.nsegments{ui.selecteddataset}(ui.selectedneuron)
            
            oerec(ui.selectedneuron).data(segmentind).info.(fieldname) = v;
            
        end
        oedb = modify_neurondata(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, oerec(ui.selectedneuron));
        
    case 'segment'
        
        oerec(ui.selectedneuron).data(ui.selectedsegment).info.(fieldname) = v;
        oedb = modify_neurondata(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, oerec(ui.selectedneuron));
        
end
setappdata(handles.oedatabrowser_figure,'oedb',oedb);


function settingsindex_edit_Callback(hObject, eventdata, handles)
setdatainfofield('settingsindex', hObject, handles);


function settingsindex_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function file_exportdataset_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, algind, algname, oedb] = oedb_getbasics(handles);
if oedb.issimulation(ui.selecteddataset)
    error('not yet implemented');
end
oerec = fetch_dataset_oerecarray(oedb, fileinfo, ui.selecteddataset); %#ok<NASGU>
[f, p] = uiputfile('*.mat','Export dataset',oedb.sourcedatafiles{ui.selecteddataset});
if isnumeric(f), return; end
save([p f],'oerec','-mat');


function focalplaneindex_edit_Callback(hObject, eventdata, handles)
setdatainfofield('focalplaneindex', hObject, handles);


function focalplaneindex_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function sbm_exporttofile_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles);
if ~strcmpi(algname, 'sbm'), return; end
[ff, pp] = uiputfile('*.gpp', 'Export GPU task problem');
if isnumeric(ff), return; end
true_st = []; %known APs not yet supported on GPU
q_spike = [];

[f,~,dt,~,~,nA2D,] = oedb_fetch_timeseries(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment); %f,it column vectors; ns row vector
results = fetch_results(oedb, fileinfo, ui.selecteddataset, n, s, algind);
if isempty(fieldnames(results.params)), return; end

sbm.toraw([pp ff], f{1}, dt, results.opts, results.params, true_st, nA2D, q_spike);


function showremovedf_Callback(hObject, eventdata, handles)
ui_toggleon(hObject);
show_data(handles);

function sbm_clearneuronspecificP_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles);
if ~ui.selecteddataset, return; end
if strcmpi(selectiontype, 'dataset')
    
    neuronlist = (1:oedb.nneurons(ui.selecteddataset))';
    [neuronlist,ok] = listdlg('ListString',oedb.neuronnames{ui.selecteddataset},'selectionmode','multiple','name','Clear neuron-specific data','InitialValue',neuronlist,'ListSize',[600 600]);
    if ~ok || isempty(neuronlist), return; end
    
else
    
    neuronlist = ui.selectedneuron;
    
end

neuron_specific_fields = {'S' 'gain' 'zeta' 'fdc' 'fr' 'R' 'mu_r_init' 'sigma_r_init'}; %FIXME get this from somewhere else

for n = reshape(neuronlist, 1, [])
    
    if ~strcmpi(selectiontype, 'segment')
        segmentlist = 1:oedb.nsegments{ui.selecteddataset}(n);
    end
    for s = segmentlist
        
        results = fetch_results(oedb, fileinfo, ui.selecteddataset, n, s, algind);
        
        for fi = 1:numel(neuron_specific_fields)
            
            if isfield(results.params, neuron_specific_fields{fi})
                
                results.params = rmfield(results.params, neuron_specific_fields{fi});
                
            end
            
        end
        
        oedb = oedb_assign_results(oedb, fileinfo, ui.selecteddataset, n, s, algind, results);
        
    end
    
end
setappdata(handles.oedatabrowser_figure,'oedb',oedb);


function sbm_fswin_explore_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles);
if ~ui.selecteddataset || ~strcmpi(algname, 'sbm'), return; end

if ~oedb.opt.storedatainmemory %fixme remove this restriction
    
    uiwait(errordlg('Data must be stored in memory to explore fs winsize', 'Error', 'modal'));
    return;
    
end
[~,it,~,~] = oedb_fetch_timeseries(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment);
dt = median(diff(it));

mi = ui_getval(handles.sbm_iter_edit);
ui_setval(handles.sbm_iter_edit, 0);

answer = inputdlg({'Number of repeats (3+)', 'Mininum window size (sec)', 'Maximum window size (sec)', 'Window size step (sec)'}, ...
    'Filter-smoother window size exploration', 1, {'25', '0.0', '1.0', num2str(dt)}); if isempty(answer), return; end
nrepeats = str2double(answer{1});
winsizes = str2double(answer{2}):str2double(answer{4}):str2double(answer{3}) - sqrt(eps); %will be rounded up
if isempty(nrepeats) || isnan(nrepeats) || imag(nrepeats) ~= 0 || nrepeats <= 2 || nrepeats ~= round(nrepeats), return; end

wb = waitbar(0, 'Exploring filter-smoother window size');
[corr, dr, fp] = deal(nan(nrepeats, numel(winsizes)));
for wind = 1:numel(winsizes)
    try
        waitbar((wind - 1) / numel(winsizes), wb); drawnow;
    catch
    end
    
    for k = 1:nrepeats
        
        oedb2 = oedb;
        
        switch selectiontype
            
            case 'dataset'
                
                for neuronind = 1:oedb.nneurons(ui.selecteddataset)
                    
                    oedb2 = oedatabrowser_detectspikes(oedb2, handles, fileinfo, algind, algname, ui.selecteddataset, neuronind); %all segments
                    
                end
                
            otherwise
                
                oedb2 = oedatabrowser_detectspikes(oedb2, handles, fileinfo, algind, algname, ui.selecteddataset, ui.selectedneuron, segmentlist);
                
        end
        
        corr(k, wind) = oedb2.stats.bysegment{ui.selecteddataset}{ui.selectedneuron}(algind, ui.selectedsegment).corr_spiketimes;
        dr(k, wind)   = oedb2.stats.bysegment{ui.selecteddataset}{ui.selectedneuron}(algind, ui.selectedsegment).det_spiketimes;
        fp(k, wind)   = oedb2.stats.bysegment{ui.selecteddataset}{ui.selectedneuron}(algind, ui.selectedsegment).fp_spiketimes;
        
    end
    
end
if isvalid(wb); close(wb); end
assignin('base','corr',corr);
assignin('base','dr',dr);
assignin('base','fp',fp);
assignin('base','winsizes',winsizes);
ui_setval(handles.sbm_iter_edit, mi);


function sbm_vartest_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles);
if ~ui.selecteddataset || ~strcmpi(algname, 'sbm'), return; end

if ~oedb.opt.storedatainmemory %fixme remove this restriction
    
    uiwait(errordlg('Data must be stored in memory to test variance', 'Error', 'modal'));
    return;
    
end

answer = inputdlg({'Number of repeats (3+)'}, 'Test Monte Carlo variance', 1, {'10'}); if isempty(answer), return; end
nrepeats = str2double(answer{1});
if isempty(nrepeats) || isnan(nrepeats) || imag(nrepeats) ~= 0 || nrepeats <= 2 || nrepeats ~= round(nrepeats), return; end
if isempty(segmentlist)
    segmentlist = 1:oedb.nsegments{ui.selecteddataset}(ui.selectedneuron);
end

mi = ui_getval(handles.sbm_iter_edit);
ui_setval(handles.sbm_iter_edit, 0);

set(handles.oedatabrowser_figure,'pointer','watch'); drawnow;
L = zeros(1, nrepeats);
wb = waitbar(0, 'Testing Monte Carlo variance');
logcondpnextobs = [];
[corr, dr, fp] = deal(nan(1, nrepeats));
ns = [];
for k = 1:nrepeats
    try
        waitbar((k - 1) / nrepeats, wb); drawnow;
    catch
    end
    
    %fixme find a way to turn off gfit when doing this
    oedb2 = oedatabrowser_detectspikes(oedb, handles, fileinfo, algind, algname, ui.selecteddataset, ui.selectedneuron, segmentlist);
    logcondpnextobs_next = [];
    
    ns_next = [];
    for j = 1:numel(segmentlist)
        
        lik = oedb2.results{algind, ui.selecteddataset}{ui.selectedneuron}(segmentlist(j)).outputvars.lik;
        L(k) = L(k) + lik.logmarglik;
        logcondpnextobs_next = [logcondpnextobs_next lik.logcondpnextobs];
        ns_next = cat(2, ns_next, oedb2.results{algind, ui.selecteddataset}{ui.selectedneuron}(segmentlist(j)).spikecounts);
        
    end
    ns = cat(1, ns, ns_next);
    
    if strcmpi(selectiontype, 'segment')
        
        corr(k) = oedb2.stats.bysegment{ui.selecteddataset}{ui.selectedneuron}(algind, ui.selectedsegment).corr_spiketimes;
        dr(k)   = oedb2.stats.bysegment{ui.selecteddataset}{ui.selectedneuron}(algind, ui.selectedsegment).det_spiketimes;
        fp(k)   = oedb2.stats.bysegment{ui.selecteddataset}{ui.selectedneuron}(algind, ui.selectedsegment).fp_spiketimes;
        
    elseif strcmpi(selectiontype, 'neuron')
        
        corr(k) = oedb2.stats.byneuron{ui.selecteddataset}(algind, ui.selectedneuron).corr_spiketimes;
        dr(k)   = oedb2.stats.byneuron{ui.selecteddataset}(algind, ui.selectedneuron).det_spiketimes;
        fp(k)   = oedb2.stats.byneuron{ui.selecteddataset}(algind, ui.selectedneuron).fp_spiketimes;
        
    elseif strcmpi(selectiontype, 'dataset')
        
        error('variance testing of full datasets is not yet supported');
        
    else
        
        error('invalid selectiontype: %s', selectiontype);
        
    end
    
    logcondpnextobs = [logcondpnextobs; logcondpnextobs_next];
    
end
assignin('base','vartest',orderfields(struct('logcondpnextobs', logcondpnextobs, 'corr', corr, 'dr', dr, 'fp', fp,'ns', ns, 'L', L, 'ui', ui, 'odbfile', fileinfo.odbfile)));
if isvalid(wb); close(wb); end
ui_setval(handles.sbm_iter_edit, mi);
set(handles.oedatabrowser_figure,'pointer','arrow');

figh = figure(); ax = [];
ax(1, end + 1) = subplot(4, 1, 1, 'parent',figh);
plot(L,'k','parent',ax(end));
hold on;
XL = get(ax(end),'xlim');
plot(XL', mean(L) * [1 1]','b','parent',ax(end));
plot(XL', mean(L) + std(L) * [-1 1; -1 1], 'g','parent',ax(end));
title(ax(end), sprintf('N %g fswin % g Mean %g Var %g SD %g',ui_getval(handles.nparticles_edit), ui_getval(handles.fswin_edit), mean(L), var(L, 1), std(L, 1)));
ylabel('log. P(F)');

ax(1, end + 1) = subplot(4, 1, 2, 'parent', figh);
plot(corr,'k','parent',ax(end));
hold on;
plot(XL', mean(corr) * [1 1]','b','parent',ax(end));
plot(XL', mean(corr) + std(corr) * [-1 1; -1 1], 'g','parent',ax(end));
ylabel('corr');

ax(1, end + 1) = subplot(4, 1, 3, 'parent', figh);
plot(dr,'k','parent',ax(end));
hold on;
plot(XL', mean(dr) * [1 1]','b','parent',ax(end));
plot(XL', mean(dr) + std(dr) * [-1 1; -1 1], 'g','parent',ax(end));
ylabel('detrate');

ax(1, end + 1) = subplot(4, 1, 4, 'parent', figh);
plot(fp,'k','parent',ax(end));
hold on;
plot(XL', mean(fp) * [1 1]','b','parent',ax(end));
plot(XL', mean(fp) + std(fp) * [-1 1; -1 1], 'g','parent',ax(end));
ylabel('fprate');
set(ax,'tickdir','out','box','off');

function sbm_profileFR_Callback(hObject, eventdata, handles)
ui_toggleon(hObject);


function sbm_profilezeta_Callback(hObject, eventdata, handles)
ui_toggleon(hObject);


function maxjitter_edit_Callback(hObject, eventdata, handles)


function maxjitter_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function roc_sdnorm_Callback(hObject, eventdata, handles)
ui_toggleon(hObject);


function roc_binarize_Callback(hObject, eventdata, handles)
ui_toggleon(hObject);


function meanfvisi_Callback(hObject, eventdata, handles)
[ui, fileinfo, ~, ~, ~, oedb] = oedb_getbasics(handles);
w = 0.5;
meanfvstdiff(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, w);


function firingratesplot_Callback(hObject, eventdata, handles)
[ui, ~, ~, ~, ~, oedb] = oedb_getbasics(handles);
oedb_firingrateplot(oedb, ui.selecteddataset);


function segmentlist = get_selected_segments(oedb, ui, selectiontype, neuron)
if strcmpi(selectiontype, 'segment')
    
    segmentlist = ui.selectedsegment;
    
else
    
    segmentlist = 1:oedb.nsegments{ui.selecteddataset}(neuron);
    
end


function spikes2bob_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles);
p = uigetdir;
if isnumeric(p)
    
    return;
    
end
if p(end) ~= filesep
    
    p = [p filesep];
    
end
neuronlist = get_selected_neurons(oedb, ui, selectiontype);
ds = ui.selecteddataset;

%fixme this is very inefficient if there are multiple ROIs in the same bobfile
[nsegtotal, nsegfound] = deal(0);
for n = neuronlist
    
    oerec = fetch_neurondata(oedb, fileinfo, ds, n);
    
    segmentlist = get_selected_segments(oedb, ui, selectiontype, n);
    nsegtotal = nsegtotal + numel(segmentlist);
    
    for s = segmentlist
        
        r = fetch_results(oedb, fileinfo, ds, n, s, algind);
        if any(isnan(r.spiketimes_window)), continue; end
        
        bobfile = [p oerec.data(s).imagefilepartialpath];
        
        if exist(bobfile, 'file') ~= 2
            
            [pp, ff, ee] = fileparts(oerec.data(s).imagefilepartialpath);
            bobfile = [p ff ee];
            if exist(bobfile, 'file') ~= 2
                
                continue; %file not found
                
            end
            
        end
        
        bobdata = load_bobdata(bobfile);
        oerec_roi_pixel_indices = sort(find(oerec.data(s).roi));
        
        avh = bobdata.image.height;
        if bobdata.align.inuse
            
            map = unique(bobdata.align.map_assoc(oerec.data.frameindex));
            if numel(map) ~= 1, continue; end
            
            if map ~= 0
                
                avh = bobdata.align.mapdims(map, 1);
                
            end
            
        end
        
        bob_roi_pixel_indices = cellfun(@(v) sort(v(:,2) + (v(:,1) - 1) * avh), bobdata.roi.interior, 'uniformoutput', false);
        
        bob_roi_index = find(cellfun(@(v) numel(v) == numel(oerec_roi_pixel_indices) && all(v == oerec_roi_pixel_indices), bob_roi_pixel_indices));
        
        if numel(bob_roi_index) ~= 1, continue; end
        
        [~,it,dt] = extract_data_fromoedatasegment(oerec.data(s));
        edges = [it - dt / 2; it(end) + dt / 2];
        
        hc = histc(r.spiketimes, edges);
        bobdata.roi.ns(:, bob_roi_index) = hc(1:end - 1);
        bobdata.roi.ns(it < r.spiketimes_window(1) | it > r.spiketimes_window(2), bob_roi_index) = nan;
        
        nsegfound = nsegfound + 1;
        
        bobdata = save_bobdata(bobdata, bobfile);
        close_bobdata(bobdata);
        
    end
    
end
fprintf('Successfully assigned spike times to bob files for %d / %d data segments across %d neurons\n', nsegfound, nsegtotal, numel(neuronlist));


function scalefactorsdnorm_Callback(hObject, eventdata, handles)
ui_toggleon(hObject);


function train_currentds_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles);
oerec = fetch_dataset_oerecarray(oedb, fileinfo, ui.selecteddataset);
ds = ui.selecteddataset;
opts = optsfromui(handles);

if ui_ison(handles.train_crossval)
    
    wb = waitbar(0, sprintf('Training %s with cross validation', algname));
    for n = 1:oedb.nneurons(ds)
        
        params = oedb.trainfunctions{algind}(oerec([1:(n - 1), n+1:end]), opts); %leave-one-out cross-validation
        oedb = assign_params_to_neuron(params, algind, ds, n, oedb, fileinfo);
        try
            waitbar(n / oedb.nneurons(ds), wb);
        catch
        end
        
    end
    if isvalid(wb); close(wb); end
    
else
    
    params = oedb.trainfunctions{algind}(oerec, opts);
    
    for n = 1:oedb.nneurons(ds)
        
        oedb = assign_params_to_neuron(params, algind, ds, n, oedb, fileinfo);
        
    end
    
end
setappdata(handles.oedatabrowser_figure,'oedb',oedb);


function oedb = assign_params_to_neuron(params, algind, dataset, neuron, oedb, fileinfo)
for s = 1:oedb.nsegments{dataset}(neuron)
    
    r = fetch_results(oedb, fileinfo, dataset, neuron, s, algind);
    r.params = params;
    oedb = oedb_assign_results(oedb, fileinfo, dataset, neuron, s, algind, r);
    
end


function train_crossval_Callback(hObject, eventdata, handles)
ui_toggleon(hObject);


function crosscorr_Callback(hObject, eventdata, handles)
[ui, fileinfo, ~, ~, ~, oedb, ~] = oedb_getbasics(handles);
datasetind = ui.selecteddataset;
resultspresent = resultspresent_eachalg(datasetind, oedb);
if isempty(resultspresent), return; end

[ii, ok] = listdlg('ListString',oedb.algnames(resultspresent),'selectionmode','multiple','name','Cross-correlation','InitialValue',1:numel(resultspresent),'ListSize',[600 600]);
if ~ok || isempty(ii), return; end

answer = inputdlg({'Shift step (ms)', 'Max shift (ms)', 'Smoothing sigma (ms)'}, 'Cross correlation', 1, {'1', '200', '50'});
if isempty(answer), return; end
shiftstep_s = str2double(answer{1}) / 1000;
maxshift_s = str2double(answer{2}) / 1000;
sigma_s = str2double(answer{3}) / 1000;

crosscorr = cell(oedb.nneurons(datasetind), numel(ii));
wb = waitbar(0, 'computing cross correlation');
for k = 1:numel(ii)
    
    algind = ii(k);
    for neuronind = 1:oedb.nneurons(datasetind)
        
        fracdone = (neuronind + (k - 1) * oedb.nneurons(datasetind)) / (numel(ii) * oedb.nneurons(datasetind));
        try
            waitbar(fracdone, wb); drawnow;
        catch
        end
        
        [crosscorr{neuronind, k}, shift_s] = oedb_calc_crosscorr(maxshift_s, shiftstep_s, sigma_s, oedb, fileinfo, algind, datasetind, neuronind);
        
    end
    
end
if isvalid(wb); close(wb); end

co = {'k' 'c' 'r' 'y' 'g' 'm' 'b'};
figure('name',sprintf('sigma = %f ms', sigma_s * 1000));
for k = 1:numel(ii)
    
    subplot(numel(ii) + 1, 1, k);
    nextcolor = co{mod(k - 1, numel(co)) + 1};
    cc = cat(2, crosscorr{:, k});
    mc = mean_nonnan(cc, 2);
    plot(-shift_s, cc, 'color', nextcolor);
    hold on;
    plot(shift_s, mc, 'color', nextcolor, 'linewidth', 2);
    ylabel(oedb.algnames{ii(k)});
    
    subplot(numel(ii) + 1, 1, numel(ii) + 1);
    plot(-shift_s, mc, 'color', nextcolor);
    hold on;
    
end
set(gca,'tickdir','out','box','off');

assignin('base', 'ccr', struct('crosscorr', {crosscorr}, 'shift_s', {shift_s}, 'sigma_s', sigma_s));


function resultspresent = resultspresent_eachalg(dataset, oedb)
resultspresent = [];
for algind = 1:numel(oedb.algnames)
    if ~isempty(oedb.stats.bydataset(algind, dataset).corr_spikecounts) && ~isnan(oedb.stats.bydataset(algind, dataset).corr_spikecounts)
        resultspresent = [resultspresent algind];
    end
end


function singlestimingerrorhist_Callback(hObject, eventdata, handles)
[ui, fileinfo, ~, ~, ~, oedb, ~] = oedb_getbasics(handles);
dataset = ui.selecteddataset;
resultspresent = resultspresent_eachalg(dataset, oedb);
if isempty(resultspresent), return; end

[ii, ok] = listdlg('ListString',oedb.algnames(resultspresent),'selectionmode','multiple','name','Timing error for single APs','InitialValue',1:numel(resultspresent),'ListSize',[600 600]);
if ~ok || isempty(ii), return; end

answer = inputdlg( ...
    {'Time before spike [ms]' 'Time after spike [ms]' 'Bin size [ms]' 'Min. number of true single APs' 'Min. number of inferred APs' 'Max. spont. firing rate (Hz)'}, ...
    'Timing error for single APs', 1, {'500' '500' '25' '5' '5' '1'});
if isempty(answer), return; end
hist_win = [-str2double(answer{1}) str2double(answer{2})] / 1000;
hist_binsize = str2double(answer{3}) / 1000;
minspikes_true = str2double(answer{4});
minspikes_inferred = str2double(answer{5});
frmax = str2double(answer{6});

unitless = ismember(oedb.algnames(ii), {'CFOOPSI', 'FOOPSI'});

[psth, nspikes_true, nspikes_inferred, ok] = deal(cell(1, numel(ii)));
[meantdiff, bias] = deal(nan(numel(ii), oedb.nneurons(dataset)));
[meantdiff_mean, meantdiff_sd, bias_mean, bias_sd] = deal(nan(1, numel(ii)));
for k = 1:numel(ii)
    
    algind = resultspresent(ii(k));
    [psth{k}, bc, nspikes_true{k}, meantdiff(k, :), ~, bias(k, :)] = oedb_estspikes_psth(oedb, fileinfo, dataset, algind, hist_binsize, hist_win);
    nspikes_inferred{k} = nspikes_true{k} .* sum(psth{k});
    
    fr = cat(2, oedb.stats.byneuron{dataset}(k, :).fr_true_spiketimes);
    
    ok{k} = nspikes_true{k} >= minspikes_true & fr <= frmax;
    if ~unitless(k)
        ok{k} = ok{k} & nspikes_inferred{k} >= minspikes_inferred;
    end
    
    y = meantdiff(k, :);
    y(~ok{k}) = [];
    y(isnan(y)) = [];
    meantdiff_mean(k) = mean(y);
    meantdiff_sd(k)   = std(y);
    
    y = bias(k, :);
    y(~ok{k}) = [];
    y(isnan(y)) = [];
    bias_mean(k) = mean(y);
    bias_sd(k)   = std(y);
    
end

assignin('base', ...
    'timingtest',orderfields(struct('psth',{psth},'bc',bc,'nspikes',{nspikes_true},'ok',{ok}, ...
    'meantdiff',meantdiff,'meantdiff_mean',meantdiff_mean,'meantdiff_sd',meantdiff_sd, ...
    'bias',meantdiff,'bias_mean',bias_mean,'bias_sd',bias_sd ...
    )));

[~, si] = sort(meantdiff_mean);
mm = nan(numel(bc), numel(ii));
for k = 1:numel(ii)
    
    if unitless(k)
        
        psth{k}(:, ok{k}) = bsxfun(@rdivide, psth{k}(:, ok{k}), sum(psth{k}(:, ok{k}), 1));  %normalize so the psth sums to 1
        
    end
    mm(:, k) = mean_nonnan(psth{k}(:, ok{k}), 2);
    
end

ff = figure;
co = {'k' 'c' 'r' 'y' 'g' 'm' 'b'};
for v = 1:numel(ii)
    
    u = si(v);
    subplot(numel(ii) + 1, 1, v, 'parent', ff);
    nextcolor = co{mod(u - 1, numel(co)) + 1};
    
    yd = psth{u}(:, ok{u});
    if unitless(u)
        
        yd = bsxfun(@rdivide, yd, max(yd, [], 1));
        
    end
    
    plot(bc, yd, nextcolor);
    
    ylabel([oedb.algnames{resultspresent(ii(u))} ' n = ' num2str(sum(ok{u}))]);
    if v < numel(ii)
        
        set(gca,'xtick',[]);
        
    end
    
end
subplot(numel(ii) + 1, 1, numel(ii) + 1,'parent',ff);
for u = 1:numel(ii)
    
    nextcolor = co{mod(u - 1, numel(co)) + 1};
    plot(bc, mm(:, u), nextcolor);  hold on;
    
end
ylabel('averages');
ax = findobj(gcf,'type','axes');
set(ax,'tickdir','out','box','off','color','none');

figure
subplot(2, 1, 1);
bar(meantdiff_mean(si));
for v = 1:numel(ii)
    
    u = si(v);
    line([v; v], meantdiff_mean(u) + [0; meantdiff_sd(u)]);
    
end
ylabel('Timing error');
subplot(2, 1, 2);
bar(bias_mean(si));
for v = 1:numel(ii)
    
    u = si(v);
    line([v; v], bias_mean(u) + [0; bias_sd(u)] * sign(bias_mean(u)));
    
end
ylabel('Bias');
ax = findobj(gcf,'type','axes');
set(ax,'xtick',1:numel(ii),'xticklabel',oedb.algnames(resultspresent(ii(si))));


function comparealgs_Callback(hObject, eventdata, handles)
[ui, fileinfo, ~, ~, ~, oedb, ~] = oedb_getbasics(handles);
XL = get(handles.f_axes, 'xlim');
oedb_compareresults(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment, XL);


function sbm_truevsim1ap_Callback(hObject, eventdata, handles)
[ui, fileinfo, selectiontype, algind, algname, oedb, segmentlist] = oedb_getbasics(handles);
if ~ui.selecteddataset, return; end

simvdata = sbm_simavspikeresponse_eachneuron(oedb, fileinfo, ui.selecteddataset, algind);
assignin('base','simvdata',simvdata);

for n = 1:oedb.nneurons(ui.selecteddataset)
    
    figure('name', [num2str(n) ' ' oedb.neuronnames{ui.selecteddataset}{n}]);
    plot(simvdata.t_data{n}, simvdata.dff_data{n},'k');  % data
    hold on
    plot(simvdata.t_sim{n}, simvdata.dff_sim{n},'r');  % model
    
end

figure('color','none');
subplot(2,1,1);
plot(simvdata.avval_data * 100, simvdata.avval_sim * 100, 'o')
maxval = max(max(simvdata.avval_data), max(simvdata.avval_sim)) * 100;
maxval_disp = 0.1 * ceil(maxval / 0.1);
hold on
plot([0 maxval_disp], [0 maxval_disp], 'r');
axis equal
axis([0 maxval_disp 0 maxval_disp]);
set(gca,'box','off','tickdir','out','color','none');
xlabel(sprintf('Data mean DFF %g to %g after AP', simvdata.avwin(1), simvdata.avwin(2)))
ylabel(sprintf('Simul. mean DFF %g to %g after AP', simvdata.avwin(1), simvdata.avwin(2)))
ok = ~isnan(simvdata.avval_data) & ~isnan(simvdata.avval_sim);
title(sprintf('r = %g', [0 1] * corrcoef(simvdata.avval_data(ok), simvdata.avval_sim(ok)) * [1 0]'));
subplot(2,1,2);
plot(simvdata.peakval_data * 100, simvdata.peakval_sim * 100, 'o')
maxval = max(max(simvdata.peakval_data), max(simvdata.peakval_sim)) * 100;
maxval_disp = 0.1 * ceil(maxval / 0.1);
hold on
plot([0 maxval_disp], [0 maxval_disp], 'r');
axis equal
axis([0 maxval_disp 0 maxval_disp]);
set(gca,'box','off','tickdir','out','color','none');
xlabel(sprintf('Data peak DFF %g to %g after AP', simvdata.avwin(1), simvdata.avwin(2)))
ylabel(sprintf('Simul. peak DFF %g to %g after AP', simvdata.avwin(1), simvdata.avwin(2)))
ok = ~isnan(simvdata.peakval_data) & ~isnan(simvdata.peakval_sim);
title(sprintf('r = %g', [0 1] * corrcoef(simvdata.peakval_data(ok), simvdata.peakval_sim(ok)) * [1 0]'));


function linearity_Callback(hObject, eventdata, handles)
ulcolor = [1 1 1] * .7; %color for unity lines etc.
patchsaturationratio = 0.3;
co = {'k' 'c' 'r' 'y' 'g' 'm' 'b'};

[ui, fileinfo, ~, ~, ~, oedb, ~] = oedb_getbasics(handles);
ds = ui.selecteddataset;

resultspresent = resultspresent_eachalg(ds, oedb); %list of algorithms with results for this dataset

resultspresent_withunits = resultspresent(~ismember(oedb.algnames(resultspresent), {'CFOOPSI', 'FOOPSI'}));

[ii, ok] = listdlg('ListString',oedb.algnames(resultspresent_withunits),'selectionmode','multiple','name','AP inference linearity analysis','InitialValue',1:numel(resultspresent_withunits));
if ~ok || isempty(ii), return; end
alglist = resultspresent_withunits(ii);

answer = inputdlg( ...
    {'Window size (s)' 'Window overlap fraction' 'Firing rate threshold (Hz)' 'Excluded time at recording start/end','Max. firing rate to plot (Hz)', 'Min. number of neurons for comparisons'}, ...
    'AP inference linearity analysis', 1, ...
    {'0.5' '0.9' '1.0' '0.5','50','3'});
if isempty(answer), return; end
winsize = str2double(answer{1});
overlapfrac = str2double(answer{2});
frthresh = str2double(answer{3});
edgebuffer = [1 1] * str2double(answer{4});
maxfrplots = str2double(answer{5});
maxspikes = floor(maxfrplots * winsize);
minneurons = str2double(answer{6});

frtrue = arrayfun(@(v) nanforempty(v.fr_true_spiketimes), oedb.stats.byneuron{ds});
lo = frtrue(alglist(1),:) < frthresh;

[rval, slope, yint] = deal(nan(oedb.nneurons(ds), numel(alglist)));
[m_lo, m_hi, s_lo, s_hi, nok_lo, nok_hi] = deal(nan(maxspikes + 1, numel(alglist))); %first row is 0 Hz
for u = 1:numel(alglist)
    
    mean_estspikes{u} = oedb_windowedspikes(oedb, fileinfo, ds, alglist(u), winsize, edgebuffer, overlapfrac);
    
    x = (0:size(mean_estspikes{u}, 1) - 1)' / winsize;
    y = mean_estspikes{u} / winsize;
    
    maxspikes_thisalg = min(maxspikes, size(y, 1) - 1);
    jj = 1:maxspikes_thisalg + 1;
    
    m_lo(jj, u) = mean_nonnan(y(jj, lo), 2);
    s_lo(jj, u) = sd_nonnan(y(jj, lo), 1, 2);
    nok_lo(jj, u) = sum(~isnan(y(jj, lo)), 2);
    m_hi(jj, u) = mean_nonnan(y(jj, ~lo), 2);
    s_hi(jj, u) = sd_nonnan(y(jj, ~lo), 1, 2);
    nok_hi(jj, u) = sum(~isnan(y(jj, ~lo)), 2);
    
    for k = 1:oedb.nneurons(ds)
        
        valok = ~isnan(y(:, k));
        if sum(valok) < 3, continue; end
        mb = [x(valok) ones(sum(valok), 1)] \ y(valok, k);
        slope(k, u) = mb(1);
        yint(k, u) = mb(2);
        cmat = corrcoef(x(valok), y(valok, k));
        rval(k, u) = cmat(2, 1);
        
    end
    
end

fig_all = figure('name', sprintf('linearity analysis, all algorithms together, window size %g s, overlap fraction %g', winsize, overlapfrac));
ax_all_lo = subplot(2,2,1,'parent', fig_all, 'nextplot','add', 'dataaspectratio', [1 1 1], 'tickdir', 'out', 'box', 'off');
plot([0 maxfrplots], [0 maxfrplots], 'color', ulcolor, 'parent', ax_all_lo);
ax_all_hi = subplot(2,2,2,'parent', fig_all, 'nextplot','add', 'dataaspectratio', [1 1 1], 'tickdir', 'out', 'box', 'off');
plot([0 maxfrplots], [0 maxfrplots], 'color', ulcolor, 'parent', ax_all_hi);
ax_ratio_lo = subplot(2, 2, 3,'parent', fig_all, 'nextplot','add', 'tickdir', 'out', 'box', 'off');
plot([0 maxfrplots], [1 1], 'color', ulcolor, 'parent', ax_ratio_lo);
ax_ratio_hi = subplot(2, 2, 4,'parent', fig_all, 'nextplot','add', 'tickdir', 'out', 'box', 'off');
plot([0 maxfrplots], [1 1], 'color', ulcolor, 'parent', ax_ratio_hi);

hav_all = [];
fig_each = figure('name', sprintf('linearity analysis, all algorithms separately, window size %g s, overlap fraction %g', winsize, overlapfrac));

x = (0:maxspikes)' / winsize;
for u = 1:numel(alglist)
    
    y = mean_estspikes{u} / winsize;
    
    maxspikes_thisalg = min(maxfrplots, size(y, 1) - 1);
    jj = 1:min(maxspikes, maxspikes_thisalg) + 1;
    
    nextcolor = co{mod(alglist(u) - 1, numel(co)) + 1};
    ax_1 = subplot(numel(alglist), 1, u, 'parent', fig_each, 'nextplot','add', 'dataaspectratio', [1 1 1], 'tickdir', 'out', 'box', 'off');
    plot([0 maxfrplots], [0 maxfrplots], 'color', ulcolor, 'parent', ax_1);
    xlabel(ax_1, 'True firing rate (Hz)'); ylabel(ax_1, 'Inferred firing rate (Hz)');
    title(ax_1, oedb.algnames{alglist(u)});
    
    if any(lo)
        plot(x(jj), y(jj, lo), 'k', 'parent', ax_1);
        plot(x, m_lo(:, u), 'k', 'linewidth', 2, 'parent', ax_1);
    end
    if any(~lo)
        plot(x(jj), y(jj, ~lo), 'r', 'parent', ax_1);
        plot(x, m_hi(:, u), 'r', 'linewidth', 2, 'parent', ax_1);
    end
    
    next_m_lo = m_lo(:, u);
    next_m_lo(nok_lo(:, u) < minneurons) = nan;
    next_m_hi = m_hi(:, u);
    next_m_hi(nok_hi(:, u) < minneurons) = nan;
    
    [hav_all(end + 1), pp_lo] = msdplot(x, next_m_lo, next_m_lo + s_lo(:, u), next_m_lo - s_lo(:, u), ax_all_lo, nextcolor);
    [~,                pp_hi] = msdplot(x, next_m_hi, next_m_hi + s_hi(:, u), next_m_hi - s_hi(:, u), ax_all_hi, nextcolor);
    cc = get(hav_all(end), 'color');
    patchcolor = max(0, min(1, 1 - (1 - cc) * patchsaturationratio));
    
    [~, pp_lo_ratio] = msdplot(x, next_m_lo ./ x, (next_m_lo + s_lo(:, u)) ./ x, (next_m_lo - s_lo(:, u)) ./ x, ax_ratio_lo, nextcolor);
    [~, pp_hi_ratio] = msdplot(x, next_m_hi ./ x, (next_m_hi + s_hi(:, u)) ./ x, (next_m_hi - s_hi(:, u)) ./ x, ax_ratio_hi, nextcolor);
    
    set([pp_lo pp_hi pp_lo_ratio pp_hi_ratio], 'facecolor', patchcolor);
    
    XL = get(ax_1,'xlim');
    set(ax_1,'xlim',[0 min(XL(2), maxfrplots)]);
    
end
ymax = max(get(ax_all_lo, 'ylim') * [0 1]', get(ax_all_hi, 'ylim') * [0 1]');
for nextax = [ax_all_lo ax_all_hi]
    XL = get(nextax,'xlim');
    set(nextax,'xlim',[0 min(XL(2), maxfrplots)], 'ylim', [0 ymax]);
end
ymax = max(get(ax_ratio_lo, 'ylim') * [0 1]', get(ax_ratio_hi, 'ylim') * [0 1]');
for nextax = [ax_ratio_lo ax_ratio_hi]
    XL = get(nextax,'xlim');
    set(nextax,'xlim',[1 / winsize min(XL(2), maxfrplots)], 'ylim', [0 ymax]);
end

legend(hav_all, oedb.algnames{alglist});
xlabel(ax_all_lo, 'True firing rate (Hz)'); ylabel(ax_all_lo, 'Inferred firing rate (Hz)');
xlabel(ax_all_hi, 'True firing rate (Hz)'); ylabel(ax_all_hi, 'Inferred firing rate (Hz)');
xlabel(ax_ratio_lo, 'True firing rate (Hz)'); ylabel(ax_ratio_lo, 'Inferred / true firing rate ratio');
xlabel(ax_ratio_hi, 'True firing rate (Hz)'); ylabel(ax_ratio_hi, 'Inferred / true firing rate ratio');

graphs_figure = figure('name',sprintf('linearity analysis, graphs, window size %g s, overlap fraction %g', winsize, overlapfrac));
for u = 1:numel(alglist)
    
    nextcolor = co{mod(alglist(u) - 1, numel(co)) + 1};
    
    subplot(2, 2, 1);
    y = slope(lo, u);
    y(isnan(y)) = [];
    if ~isempty(y)
        bar(u, mean(y), .5, 'facecolor', nextcolor, 'edgecolor', 'k'); hold on;
        plot([1; 1] * u, mean(y) + [0; std(y, 1)], 'k');
        plot(u + .4, y, 'o', 'markerfacecolor', nextcolor, 'markeredgecolor', 'k');
    end
    ylabel('Slope of inferred vs. true firing rates');
    title('putative pyramidals');
    
    subplot(2, 2, 2);
    y = slope(~lo, u);
    y(isnan(y)) = [];
    if ~isempty(y)
        bar(u, mean(y), .5, 'facecolor', nextcolor, 'edgecolor', 'k'); hold on;
        plot([1; 1] * u, mean(y) + [0; std(y, 1)], 'k');
        plot(u + .4, y, 'o', 'markerfacecolor', nextcolor, 'markeredgecolor', 'k');
    end
    ylabel('Slope of inferred vs. true firing rates');
    title('putative interneurons');
    
    subplot(2, 2, 3);
    y = rval(lo, u) .^ 2;
    y(isnan(y)) = [];
    if ~isempty(y)
        bar(u, mean(y), .5, 'facecolor', nextcolor, 'edgecolor', 'k'); hold on;
        plot([1; 1] * u, mean(y) + [0; std(y, 1)], 'k');
        plot(u + .4, y, 'o', 'markerfacecolor', nextcolor, 'markeredgecolor', 'k');
    end
    ylabel('Regression R^2');
    title('putative pyramidals');
    
    subplot(2, 2, 4);
    y = rval(~lo, u) .^ 2;
    y(isnan(y)) = [];
    if ~isempty(y)
        bar(u, mean(y), .5, 'facecolor', nextcolor, 'edgecolor', 'k'); hold on;
        plot([1; 1] * u, mean(y) + [0; std(y, 1)], 'k');
        plot(u + .4, y, 'o', 'markerfacecolor', nextcolor, 'markeredgecolor', 'k');
    end
    ylabel('Regression R^2');
    title('putative interneurons');
    
end
set(findobj(graphs_figure,'type','axes'), 'xtick', 1:numel(alglist), 'xticklabel', oedb.algnames(alglist));

assignvars = {'m_lo' 'm_hi' 's_lo' 's_hi' 'mean_estspikes' 'slope' 'yint' 'rval' 'edgebuffer' 'frthresh' 'nok_lo' 'nok_hi'};
lintest = struct('algnames', {oedb.algnames(alglist)}, 'dataset', oedb.datasetnames{ds}, 'odbfile', fileinfo.odbfile);
for j = 1:numel(assignvars);
    
    lintest.(assignvars{j}) = eval(assignvars{j});
    
end
assignin('base','lintest',lintest);


function burstinf_vs_isi_Callback(hObject, eventdata, handles)
ulcolor = [1 1 1] * .7; %color for unity lines etc.
patchsaturationratio = 0.3;
co_algs = {'k' 'c' 'r' 'y' 'g' 'm' 'b'};
co_nspikes = {'k' 'b' 'g' 'r' 'c' 'm' 'y'};
binsize_isihist = 1e-3;

[ui, fileinfo, ~, ~, ~, oedb, ~] = oedb_getbasics(handles);
ds = ui.selecteddataset;

resultspresent = resultspresent_eachalg(ds, oedb); %list of algorithms with results for this dataset

resultspresent_withunits = resultspresent(~ismember(oedb.algnames(resultspresent), {'CFOOPSI', 'FOOPSI'}));
if isempty(resultspresent_withunits), return; end

[ii, ok] = listdlg('ListString',oedb.algnames(resultspresent_withunits),'selectionmode','multiple','name','Burst inference accuracy vs. ISI','InitialValue',1:numel(resultspresent_withunits));
if ~ok || isempty(ii), return; end
alglist = resultspresent_withunits(ii);

answer = inputdlg( ...
    {'Max. ISI within burst (ms)' 'Min. time without APs before burst (ms)' 'Min. time without APs after burst (ms)' 'Max. time difference to inferred APs (ms)' 'Max. APs in burst' 'Mean ISI bin size (ms)' 'Firing rate threshold (Hz)'}, ...
    'Burst inference accuracy vs. ISI', 1, ...
    {'100' '400' '400' '200','3', '10', '1.0'});
if isempty(answer), return; end
maxisi = str2double(answer{1}) / 1000;
mintbefore = str2double(answer{2}) / 1000;
mintafter = str2double(answer{3}) / 1000;
maxtdiff = str2double(answer{4}) / 1000;
maxaps = str2double(answer{5});
dISI = str2double(answer{6}) / 1000;
frthresh = str2double(answer{7});

frtrue = arrayfun(@(v) nanforempty(v.fr_true_spiketimes), oedb.stats.byneuron{ds});
lo = frtrue(alglist(1),:) < frthresh;

ed = dISI * (0:ceil(maxisi / dISI));
bc = ed(1:end - 1) + diff(ed(1:2)) / 2;
ed_isihist = binsize_isihist * (0:ceil(maxisi / binsize_isihist));
bc_isihist = ed_isihist(1:end - 1) + diff(ed_isihist(1:2)) / 2;

true_est_meanisi = cell(1, numel(alglist));
meaninf_each = nan(numel(bc), maxaps, numel(alglist), oedb.nneurons(ds));
[meaninf_av_lo, meaninf_sd_lo, meaninf_av_hi, meaninf_sd_hi] = deal(nan(numel(bc), maxaps, numel(alglist)));
for u = 1:numel(alglist)
    
    %first column is number of true APs, second is number of inferred APs, third is mean ISI
    true_est_meanisi{u} = oedb_burstinf_vs_n_and_meanisi(oedb, fileinfo, ds, alglist(u), maxisi, mintbefore, mintafter, maxtdiff, maxaps);
    
    for naps = 2:maxaps
        for ibin = 1:numel(bc)
            
            for ineuron = 1:oedb.nneurons(ds)
                
                ii = true_est_meanisi{u}{ineuron}(:, 3) >= ed(ibin) & true_est_meanisi{u}{ineuron}(:, 3) < ed(ibin + 1) & true_est_meanisi{u}{ineuron}(:, 1) == naps;
                meaninf_each(ibin, naps, u, ineuron) = mean(true_est_meanisi{u}{ineuron}(ii, 2));
                
            end
            meaninf_av_lo(ibin, naps, u) = mean_nonnan(meaninf_each(ibin, naps, u, lo), 4);
            meaninf_sd_lo(ibin, naps, u) = sd_nonnan(meaninf_each(ibin, naps, u, lo), 1, 4);
            meaninf_av_hi(ibin, naps, u) = mean_nonnan(meaninf_each(ibin, naps, u, ~lo), 4);
            meaninf_sd_hi(ibin, naps, u) = sd_nonnan(meaninf_each(ibin, naps, u, ~lo), 1, 4);
            
        end
    end
    
end

%fixme check assumption all algs have the same data processed!
allbursts_lo = cat(1, true_est_meanisi{1}{lo});
allbursts_hi = cat(1, true_est_meanisi{1}{~lo});
isihist_lo = histc(allbursts_lo(:, 3), ed_isihist);
isihist_hi = histc(allbursts_hi(:, 3), ed_isihist);

maxinf = max(max(allbursts_lo(:, 2)), max(allbursts_hi(:, 2)));
[prob_inf_given_true_lo, prob_inf_given_true_hi] = deal(zeros(maxaps, maxinf + 1, numel(alglist))); %row is true APs, column is inferred APs + 1 (first column is for 0 APs inferred)
for u = 1:numel(alglist)
    
    for naps = 2:maxaps
        
        [nused_lo, nused_hi] = deal(0);
        for ineuron = 1:oedb.nneurons(ds)
            
            ii = true_est_meanisi{u}{ineuron}(:, 1) == naps;
            if ~any(ii), continue; end
            if lo(ineuron)
                prob_inf_given_true_lo(naps, :, u) = prob_inf_given_true_lo(naps, :, u) + reshape(histc(true_est_meanisi{u}{ineuron}(ii, 2), 0:maxinf) / sum(ii), 1, []);
                nused_lo = nused_lo + 1;
            else
                prob_inf_given_true_hi(naps, :, u) = prob_inf_given_true_hi(naps, :, u) + reshape(histc(true_est_meanisi{u}{ineuron}(ii, 2), 0:maxinf) / sum(ii), 1, []);
                nused_hi = nused_hi + 1;
            end
            
        end
        prob_inf_given_true_lo(naps, :, u) = prob_inf_given_true_lo(naps, :, u) / nused_lo;
        prob_inf_given_true_hi(naps, :, u) = prob_inf_given_true_hi(naps, :, u) / nused_hi;
        
    end
    
end

for naps = 2:maxaps
    
    isihist_lo_thiscount = histc(allbursts_lo(allbursts_lo(:,1) == naps, 3), ed_isihist);
    isihist_hi_thiscount = histc(allbursts_hi(allbursts_hi(:,1) == naps, 3), ed_isihist);
    
    fig_comp = figure('name', sprintf('ISI vs. burst inference, comparsions, %d APS, max. ISI %f', naps, maxisi));
    subplot(2,2,1,'parent', fig_comp, 'nextplot','add', 'tickdir', 'out', 'box', 'off');
    bar(bc_isihist, isihist_lo_thiscount(1:end-1));
    xlabel('Mean ISI within burst');
    ylabel('Count');
    title('Putative pyramidals');
    subplot(2,2,2,'parent', fig_comp, 'nextplot','add', 'tickdir', 'out', 'box', 'off');
    bar(bc_isihist, isihist_hi_thiscount(1:end-1));
    xlabel('Mean ISI within burst');
    ylabel('Count');
    title('Putative interneurons');
    
    ax_lo = subplot(2,2,3,'parent', fig_comp, 'nextplot','add', 'tickdir', 'out', 'box', 'off');
    xlabel('Mean ISI within burst');
    ylabel('Detected APs');
    ax_hi = subplot(2,2,4,'parent', fig_comp, 'nextplot','add', 'tickdir', 'out', 'box', 'off');
    xlabel('Mean ISI within burst');
    ylabel('Detected APs');
    hline_lo = nan(1, numel(alglist));
    for u = 1:numel(alglist)
        
        nextcolor = co_algs{mod(alglist(u) - 1, numel(co_algs)) + 1};
        
        [hline_lo(u), hpatch_lo] = msdplot(bc, meaninf_av_lo(:, naps, u), ...
            meaninf_av_lo(:, naps, u) + meaninf_sd_lo(:, naps, u), meaninf_av_lo(:, naps, u) - meaninf_sd_lo(:, naps, u), ax_lo, nextcolor);
        cc = get(hline_lo(u), 'color');
        patchcolor = max(0, min(1, 1 - (1 - cc) * patchsaturationratio));
        
        [hline_hi, hpatch_hi] = msdplot(bc, meaninf_av_hi(:, naps, u), ...
            meaninf_av_hi(:, naps, u) + meaninf_sd_hi(:, naps, u), meaninf_av_hi(:, naps, u) - meaninf_sd_hi(:, naps, u), ax_hi, nextcolor);
        
        set([hline_hi hline_lo(u)], 'linewidth', 2);        
        set([hpatch_hi hpatch_lo], 'facecolor', patchcolor);        
        
    end
    legend(hline_lo, oedb.algnames{alglist});
    
end
    
for u = 1:numel(alglist)
    
    fig_indiv = figure('name', sprintf('ISI vs. burst inference, %s, 2 to %d APS, max. ISI %f', oedb.algnames{alglist(u)}, maxaps, maxisi));
    
    subplot(3,2,1,'parent', fig_indiv, 'nextplot','add', 'tickdir', 'out', 'box', 'off');
    bar(bc_isihist, isihist_lo(1:end-1));
    xlabel('Mean ISI within burst');
    ylabel('Count');
    title('Putative pyramidals');
    
    subplot(3,2,2,'parent', fig_indiv, 'nextplot','add', 'tickdir', 'out', 'box', 'off');
    bar(bc_isihist, isihist_hi(1:end-1));
    xlabel('Mean ISI within burst');
    ylabel('Count');
    title('Putative interneurons');
    
    ax_lo = subplot(3,2,3,'parent', fig_indiv, 'nextplot','add', 'tickdir', 'out', 'box', 'off');
    xlabel('Mean ISI within burst');
    ylabel('Detected APs');
    ylim([-.5 maxinf + .5]);
    
    ax_hi = subplot(3,2,5,'parent', fig_indiv, 'nextplot','add', 'tickdir', 'out', 'box', 'off');
    xlabel('Mean ISI within burst');
    ylabel('Detected APs');
    ylim([-.5 maxinf + .5]);
    
    subplot(3,2,4,'parent', fig_indiv, 'nextplot','add', 'tickdir', 'out', 'box', 'off');
    hbar_hi = barh((0:maxinf)', prob_inf_given_true_hi(2:end, :, u)');
    ylim([-.5 maxinf + .5]);
    xlabel('Probability');
    ylabel('Inferred APs');
    title('Putative interneurons');
    
    subplot(3,2,6,'parent', fig_indiv, 'nextplot','add', 'tickdir', 'out', 'box', 'off');
    hbar_lo = barh((0:maxinf)', prob_inf_given_true_lo(2:end, :, u)');
    ylim([-.5 maxinf + .5]);
    xlabel('Probability');
    ylabel('Inferred APs');
    title('Putative pyramidals');
    
    hline_lo = nan(1, maxaps - 1);
    tsi_thisalg_lo = cat(1, true_est_meanisi{u}{lo});
    tsi_thisalg_hi = cat(1, true_est_meanisi{u}{~lo});
    
    for naps = 2:maxaps
        
        nextcolor = co_nspikes{mod(naps - 1, numel(co_nspikes)) + 1};
        set([hbar_lo(naps - 1) hbar_hi(naps - 1)], 'facecolor', nextcolor);
        
        [hline_lo(naps - 1), hpatch_lo] = msdplot(bc, meaninf_av_lo(:, naps, u), ...
            meaninf_av_lo(:, naps, u) + meaninf_sd_lo(:, naps, u), meaninf_av_lo(:, naps, u) - meaninf_sd_lo(:, naps, u), ax_lo, nextcolor);
        cc = get(hline_lo(naps - 1), 'color');
        patchcolor = max(0, min(1, 1 - (1 - cc) * patchsaturationratio));
        
        [hline_hi, hpatch_hi] = msdplot(bc, meaninf_av_hi(:, naps, u), ...
            meaninf_av_hi(:, naps, u) + meaninf_sd_hi(:, naps, u), meaninf_av_hi(:, naps, u) - meaninf_sd_hi(:, naps, u), ax_hi, nextcolor);

        ii = tsi_thisalg_lo(:, 1) == naps;
        plot(tsi_thisalg_lo(ii, 3), tsi_thisalg_lo(ii, 2), 'o', 'color', nextcolor, 'parent', ax_lo);
        ii = tsi_thisalg_hi(:, 1) == naps;
        plot(tsi_thisalg_hi(ii, 3), tsi_thisalg_hi(ii, 2), 'o', 'color', nextcolor, 'parent', ax_hi);
        
        set([hline_hi hline_lo(naps - 1)], 'linewidth', 2);        
        set([hpatch_hi hpatch_lo], 'facecolor', patchcolor);
        
    end
    legend(hline_lo, cellstr(num2str((2:maxaps)')));
       
end


function sbm_shotnoise_Callback(hObject, eventdata, handles)
ui_toggleon(hObject);
