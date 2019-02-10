function [oedb, canceled] = oedb_nonlinfit(oedb, ui, fileinfo, selectiontype, algind, handles, target_substepsize, maxtotalt_hours)
canceled = true;
%FIXME this function is a mess
savefigs = true;

parmode = (exist('matlabpool','file') && matlabpool('size') > 0) || (exist('parpool','file') && ~isempty(gcp('nocreate')));

modelname = []; %use default

if ~exist('target_substepsize', 'var') || isempty(target_substepsize) || ~exist('maxtotalt_hours', 'var') || isempty(maxtotalt_hours)
    
    ssms = nan;
    maxtotalt_hours = nan;
    while ~isfinite(ssms) || imag(ssms) ~= 0 || ssms <= 1e-6 || ~isfinite(maxtotalt_hours) || maxtotalt_hours <= 0 || imag(maxtotalt_hours) ~= 0
        
        answer = inputdlg({'Max. step size (ms)' 'Max fitting time (hours)'}, 'Nonlinear fit', 1, {'10' '40'});
        if isempty(answer), return; end
        ssms = str2double(answer{1});
        maxtotalt_hours = str2double(answer{2});
        
    end
    target_substepsize = 1e-3 * ssms;
    
end

[f,it,st,indicatorstring] = deal({});
%note that neuronind starts from one and goes up to the number of neurons used, even when some neurons are excluded. list of neuron indices in oedb are stored in neuronlist
switch selectiontype
    case 'dataset'
        
        neuronlist = (1:oedb.nneurons(ui.selecteddataset))';
        [neuronlist,ok] = listdlg('ListString',oedb.neuronnames{ui.selecteddataset},'selectionmode','multiple','name','Nonlinear fit','InitialValue',neuronlist,'ListSize',[600 600]);
        if ~ok || isempty(neuronlist), return; end
        neuronnames = oedb.neuronnames{ui.selecteddataset}(neuronlist);
        [neuronind, sessionindex, focalplaneindex] = deal([]);
        
        oerec = fetch_dataset_oerecarray(oedb, fileinfo, ui.selecteddataset);
        
        for k = reshape(neuronlist,1,[])
            
            [~, next_sessionindex, ~, ~, next_focalplaneindex] = oerec_datagroupindices_singleneuron(oerec(k));
            next_focalplaneindex(next_focalplaneindex > 0) = next_focalplaneindex(next_focalplaneindex > 0) + max([0; focalplaneindex]);
            
            focalplaneindex       = [focalplaneindex;       next_focalplaneindex];
            sessionindex = [sessionindex; next_sessionindex + max([0; sessionindex])];
            neuronind = [neuronind; max([0; neuronind]) + ones(oedb.nsegments{ui.selecteddataset}(k),1)];
            for s = 1:oedb.nsegments{ui.selecteddataset}(k)
                
                [f{end + 1, 1},it{end + 1, 1},~,~,indicatorstring{end + 1, 1},~,st{end + 1, 1},~,~,stwin] = oedb_fetch_timeseries(oedb, fileinfo, ui.selecteddataset, k, s); %f,it column vectors; ns row vector
                assert(~any(isnan(stwin)), 'Need within-frame spike times');
                
            end
            
        end
        
    case 'neuron'
        
        neuronlist = ui.selectedneuron;
        neuronind = ones(oedb.nsegments{ui.selecteddataset}(ui.selectedneuron), 1);
        neuronnames = oedb.neuronnames{ui.selecteddataset}(ui.selectedneuron);
        
        oerec = fetch_dataset_oerecarray(oedb, fileinfo, ui.selecteddataset);
        
        [~, sessionindex, ~, ~, focalplaneindex] = oerec_datagroupindices_singleneuron(oerec(ui.selectedneuron));
        for s = 1:oedb.nsegments{ui.selecteddataset}(ui.selectedneuron)
            
            [f{end + 1, 1},it{end + 1, 1},~,~,indicatorstring{end + 1, 1},~,st{end + 1, 1},~,~,stwin] = oedb_fetch_timeseries(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, s); %f,it column vectors; ns row vector
            assert(~any(isnan(stwin)), 'Need within-frame spike times');
            
        end
        
    case 'segment'
        
        neuronlist = ui.selectedneuron;
        neuronind = 1;
        neuronnames = oedb.neuronnames{ui.selecteddataset}(ui.selectedneuron);
        sessionindex = 1;
        focalplaneindex = 1;
        
        [f{1},it{1},~,~,indicatorstring{end + 1, 1},~,st{1},~,~,stwin] = oedb_fetch_timeseries(oedb, fileinfo, ui.selecteddataset, ui.selectedneuron, ui.selectedsegment); %f,it column vectors; ns row vector
        assert(~any(isnan(stwin)), 'Need within-frame spike times');
        
end
nneurons = numel(neuronlist);

if numel(f) == 1 && ui_ison(handles.onscreenonly)
    
    XL = get(handles.f_axes,'xlim') + [-1 1] * median(diff(it{1}));
    fmask = it{1} >= XL(1) & it{1} <= XL(2);
    it{1} = it{1}(fmask); f{1} = f{1}(fmask);
    
end

if any(cellfun(@isempty,indicatorstring)) || any(~strcmp(indicatorstring{1}, indicatorstring))
    
    uiwait(errordlg('Indicator must be known and the same for all fitted data', 'Error', 'modal')); return;
    
end
assert(all(strcmpi(indicatorstring{1}, indicatorstring{1})), 'all data must be from the same calcium indicator');
indicatorstring = indicatorstring{1};

basefigsavedir = [fileparts(mfilename('fullpath')) filesep '..' filesep '..' filesep 'figs' filesep];
if ~exist(basefigsavedir, 'dir') && savefigs
    mkdir(basefigsavedir);
    if ~exist(basefigsavedir, 'dir')
        
        uiwait(errordlg('Failed to create figures directory', 'Error', 'modal')); return;
        
    end
end

answer = inputdlg({'Model fit directory name'}, 'Nonlinear fit', 1, {datestr(now, 'ddmmyy')});
if isempty(answer), return; end
basefigsavedir = [basefigsavedir answer{1} filesep];
if ~exist(basefigsavedir, 'dir') && savefigs
    mkdir(basefigsavedir);
    if ~exist(basefigsavedir, 'dir')
        
        uiwait(errordlg('Failed to create figures directory', 'Error', 'modal')); return;
        
    end
end

fitmode = '';
if strcmpi(selectiontype, 'dataset')
    
    answer = questdlg('Cross-validate?', 'Model fit', 'Yes', 'No', 'Cancel', 'Yes');
    if isempty(answer) || strcmpi(answer, 'Cancel'), return; end
    if strcmpi(answer, 'Yes')
        
        fitmode = 'crossval';
        
    end
    
end
if isempty(fitmode)
    
    
    answer = questdlg('Choose fit type', 'Model fit', 'Simple', 'Multi-start from submodels', 'Multistart from file', 'Simple');
    if isempty(answer), return; end %user closed dialog without clicking a button
    if strcmpi(answer, 'Multi-start from submodels')
        
        fitmode = 'submodels';
        
    elseif strcmpi(answer, 'Multistart from file')
        
        fitmode = 'file';
        
        [ff, pp] = getfile_fromp(fileinfo.odbdir, '*.mat', 'Choose parameter file');
        if isnumeric(ff), return; end
        paramsfile = [pp ff];
        P_fromfile = load(paramsfile, 'P', '-mat');
        P_fromfile = reshape(P_fromfile.P, 1, []);
        assert(isa(P_fromfile, 'struct') && ~isempty(fieldnames(P_fromfile)) && ~isempty(P_fromfile), 'invalid file');
        
    end
    
end
always_fixed_params = {}; P_empty = [];
savefigs = false;
if strcmpi(fitmode, 'crossval')
    
    [training_neurons, testing_neurons] = crossval_groups(nneurons);
    if any(isnan(training_neurons(:))), return; end
    one_test_per_neuron = size(testing_neurons, 1) == nneurons && size(testing_neurons, 2) == 1 && numel(testing_neurons) == numel(unique(testing_neurons)) && all(unique(testing_neurons(:)) == (1:nneurons)');
    ngroups = size(training_neurons, 1);
    
    opts = struct();
    Pout_he = cell(1, ngroups);
    
    basefigsavedir_eachgroup = arrayfun(@(g) [basefigsavedir 'cv training group ' num2str(g) filesep], 1:ngroups, 'uniformoutput', false)';
    overwrite = true;
    if any(cellfun(@(v) exist(v, 'dir'), basefigsavedir_eachgroup))
        
        answer = questdlg('Target directory(s) exist. Overwrite existing files?', 'Model fit', 'Yes', 'No', 'Cancel', 'Yes');
        if isempty(answer) || strcmpi(answer, 'Cancel'), return; end
        overwrite = strcmpi(answer, 'Yes');
        
    end
    
    if parmode
        opts.verbose = 0;
        parfor g = 1:ngroups
            
            Pout_he{g} = run_crossval_group(neuronind, sessionindex, focalplaneindex, training_neurons(g, :), testing_neurons(g, :), opts, basefigsavedir_eachgroup{g}, ...
                f, it, st, neuronnames, always_fixed_params, modelname, target_substepsize, maxtotalt_hours, savefigs, overwrite);
            
        end
    else
        for g = 1:ngroups
            
            fprintf('\ntraining group %d / %d\n', g, ngroups);
            
            Pout_he{g} = run_crossval_group(neuronind, sessionindex, focalplaneindex, training_neurons(g, :), testing_neurons(g, :), opts, basefigsavedir_eachgroup{g}, ...
                f, it, st, neuronnames, always_fixed_params, modelname, target_substepsize, maxtotalt_hours, savefigs, overwrite);
            
        end
    end
    for g = 1:ngroups %assing parameters to held-out neurons
        if one_test_per_neuron
            
            k = neuronlist(testing_neurons(g));  %convert from index among used neurons (g) to index of all neurons in the selected oedb dataset (k)
            
            for s = 1:oedb.nsegments{ui.selecteddataset}(k)
                
                oedb = assign_params_toseg(oedb, fileinfo, ui.selecteddataset, k, s, algind, Pout_he{g}(s), indicatorstring);
                
            end
        end
    end
    
else
    
    if strcmpi(fitmode, 'file')
        
        opts = sbm.fit.assignopts(struct());
        params_from_file = fieldnames(P_fromfile);
        fixed_params = [always_fixed_params params_from_file(:)'];
        if ismember('k50', params_from_file)
            
            [opts.mink50_first, opts.maxk50_first, opts.min_dk50, opts.max_dk50] = deal(0, inf, 0, inf);
            
        end
        if ismember('tau_mm', params_from_file)
            
            [opts.min_tau_mm, opts.max_tau_mm] = deal(0, inf);
            
        end
        for k = 1:numel(P_fromfile)
            
            sbm.fit.globalfit( ...
                f, it, st, neuronind, neuronnames, basefigsavedir, P_fromfile(k), ...
                sessionindex, fixed_params, focalplaneindex, modelname, target_substepsize, maxtotalt_hours, savefigs, sprintf('fit from file %d.mat', k), opts);
            
        end
        
    elseif strcmpi(fitmode, 'submodels')
        
        defaultopts = sbm.fit.assignopts(struct());
        opts = struct('usek50', true, 'nbindingsteps', 4);
        P0 = struct();
        k50_first = 0.1; %uM
        maxk50_init = 0.6;
        dk50_absentstates_init = 0.0025;  %uM
        maxdk50_absentstates = 0.005;  %uM
        maxdk50_presentstates = max(max(defaultopts.max_dk50), 1.1 * (maxk50_init - k50_first));
        nsubmodels = 2^(opts.nbindingsteps - 1);
        
        submodels = dec2base(0:nsubmodels -1 , 2) == '1';
        
        if parmode
            
            opts.verbose = 0;
            parfor submodelind = 1:nsubmodels
                
                fit_from_submodel(submodels, submodelind, opts, P0, k50_first, maxk50_init, dk50_absentstates_init, maxdk50_absentstates, maxdk50_presentstates, ...
                    f, it, st, neuronind, neuronnames, basefigsavedir, sessionindex, always_fixed_params, focalplaneindex, modelname, target_substepsize, maxtotalt_hours, savefigs)
                
            end
            
        else
            
            opts.verbose = 10;
            for submodelind = 1:nsubmodels
                
                fit_from_submodel(submodels, submodelind, opts, P0, k50_first, maxk50_init, dk50_absentstates_init, maxdk50_absentstates, maxdk50_presentstates, ...
                    f, it, st, neuronind, neuronnames, basefigsavedir, sessionindex, always_fixed_params, focalplaneindex, modelname, target_substepsize, maxtotalt_hours, savefigs)
                
            end
            
        end
        %fixme assign the best one
        
    else  % simple
        
        [Pout, ~, firingrates_eachneuron] = ...
            sbm.fit.globalfit( ...
            f, it, st, neuronind, neuronnames, basefigsavedir, P_empty, ...
            sessionindex, always_fixed_params, focalplaneindex, modelname, target_substepsize, maxtotalt_hours, savefigs, 'vars.mat');
        
        Pout = sbm.fit.fit_hyperparams(Pout, neuronind, firingrates_eachneuron);
        
        for j = 1:numel(neuronlist)
            
            k = neuronlist(j);
            ii = find(neuronind == j);
            
            if strcmpi(selectiontype, 'segment')
                
                slist = ui.selectedsegment;
                
            else
                
                slist = 1:oedb.nsegments{ui.selecteddataset}(k);
                
            end
            
            for v = 1:numel(slist)
                
                s = slist(v);
                Pout(ii(v)).fr = firingrates_eachneuron(j);
                
                oedb = assign_params_toseg(oedb, fileinfo, ui.selecteddataset, k, s, algind, Pout(ii(v)), indicatorstring);
                
            end
        end
        
    end
    
end
canceled = false;


function Pout_he = run_crossval_group(neuronind, sessionindex, focalplaneindex, training_neurons, testing_neurons, opts, basefigsavedir_cv, f, it, st, neuronnames, ...
    always_fixed_params, modelname, target_substepsize, maxtotalt_hours, savefigs, overwrite)
fitmatname = 'training data fit.mat';
P_empty = [];

training_segs = ismember(neuronind, training_neurons);
[~, ~, neuronind_trainingneurons_from1] = unique(neuronind(training_segs));
training_sessions = unique(sessionindex(training_segs));
training_sessionindex_from1 = nan(sum(training_segs), 1);
for q = 1:numel(training_sessions)
    
    training_sessionindex_from1(sessionindex(training_segs) == training_sessions(q)) = q;
    
end

if exist(basefigsavedir_cv, 'dir') ~= 7
    mkdir(basefigsavedir_cv);
    assert(exist(basefigsavedir_cv, 'dir') == 7, 'failed to create directory for figures');
end
%fit all training neurons together

if overwrite || exist([basefigsavedir_cv fitmatname], 'file') ~= 2
    
    [Pout, opts, firingrates] = ...
        sbm.fit.globalfit( ...
        f(training_segs), it(training_segs), st(training_segs), neuronind_trainingneurons_from1, ...
        neuronnames(training_neurons), basefigsavedir_cv, P_empty, ...
        training_sessionindex_from1, always_fixed_params, focalplaneindex(training_segs), modelname, target_substepsize, maxtotalt_hours, savefigs, fitmatname, opts);
    
else %use existing file
    
    load([basefigsavedir_cv fitmatname], 'Pout', 'opts', 'firingrates');
    
end
%fit hyperparameters
Pout = sbm.fit.fit_hyperparams(Pout, neuronind_trainingneurons_from1, firingrates);
Pout = Pout(1);
%assign prior mode for R, S
Pout.S = exp(Pout.normal_meancov_logSR(1, 1));
Pout.R = exp(Pout.normal_meancov_logSR(2, 1));
%assign prior mean for FR (mode is zero for k <= 1)
Pout.fr = prod(Pout.gamma_kt_FR);

%fit each testing neuron individually
fixed_params_testing = union(always_fixed_params, setdiff(fieldnames(Pout), opts.neuronspecific));
for k = testing_neurons
    
    allsegs_thisneuron = neuronind == k;
    sessions_thisneuron = unique(sessionindex(allsegs_thisneuron));
    sessionindex_thisneuron_from1 = nan(sum(allsegs_thisneuron), 1);
    for q = 1:numel(sessions_thisneuron)
        
        sessionindex_thisneuron_from1(sessionindex(allsegs_thisneuron) == sessions_thisneuron(q)) = q;  % renumber sessions
        
    end
    
    varsfile = sprintf('testing data fit neuron %d.mat', k);
    
    if overwrite || exist([basefigsavedir_cv varsfile], 'file') ~= 2
        
        [Pout_he, ~, ~] = sbm.fit.globalfit( ...
            f(allsegs_thisneuron), it(allsegs_thisneuron), st(allsegs_thisneuron), ones(sum(neuronind == k), 1), ...
            neuronnames(k), basefigsavedir_cv, Pout, ...
            sessionindex_thisneuron_from1, fixed_params_testing, focalplaneindex(allsegs_thisneuron), modelname, target_substepsize, maxtotalt_hours, savefigs, varsfile, opts);
        
    else
        
        load([basefigsavedir_cv varsfile], 'Pout');
        Pout_he = Pout;
        
    end
    
end


function fit_from_submodel(submodels, submodelind, opts, P0, k50_first, maxk50_init, dk50_absentstates_init, maxdk50_absentstates, maxdk50_presentstates, ...
    f, it, st, neuronind, neuronnames, basefigsavedir, sessionindex, always_fixed_params, focalplaneindex, modelname, target_substepsize, maxtotalt_hours, savefigs)
nsubmodels = 2^(opts.nbindingsteps - 1);
if opts.verbose > 5
    
    fprintf('\nset of stable intermediates %d / %d\n', submodelind, nsubmodels);
    
end

intermediates_present = submodels(submodelind, :);
opts.max_dk50 = maxdk50_absentstates * ones(1, opts.nbindingsteps - 1); %uM
opts.max_dk50(intermediates_present) = maxdk50_presentstates;

if any(intermediates_present)
    
    P0.k50_first = k50_first;
    
else  % direct transition from apo to saturated state
    
    P0.k50_first = (k50_first + maxk50_init) / 2;
    
end

dk50 = nan(1, opts.nbindingsteps - 1);
dk50(~intermediates_present) = dk50_absentstates_init;
remaining_range = (maxk50_init - P0.k50_first) - sum(~intermediates_present) * dk50_absentstates_init;
dk50(intermediates_present) = remaining_range / sum(intermediates_present);
P0.dk50 = dk50;

varsfile = sprintf('submodel %d.mat', submodelind);

[Pfit_submodel, ~, ~] = ...
    sbm.fit.globalfit( ...
    f, it, st, neuronind, neuronnames, basefigsavedir, P0, ...
    sessionindex, always_fixed_params, focalplaneindex, modelname, target_substepsize, maxtotalt_hours, savefigs, varsfile, opts);

varsfile2 = sprintf('init from submodel %d.mat', submodelind);
if all(intermediates_present)  % no constraints to release
    
    copyfile([basefigsavedir varsfile], [basefigsavedir varsfile2], 'f');
    
else
    
    opts.max_dk50(~intermediates_present) = maxdk50_presentstates; %release the constraint and keep going
    sbm.fit.globalfit( ...
        f, it, st, neuronind, neuronnames, basefigsavedir, Pfit_submodel, ...
        sessionindex, always_fixed_params, focalplaneindex, modelname, target_substepsize, maxtotalt_hours, savefigs, varsfile2, opts);
    
end


function [training_neurons, testing_neurons] = crossval_groups(nneurons)
[training_neurons, testing_neurons] = deal(nan);
answer = {' '}; ntrain = nan;
while ~strcmpi(answer{1}, 'all') && (~isfinite(ntrain) || ntrain ~= round(ntrain) || imag(ntrain) ~= 0 || ntrain <= 0 || ntrain >= nneurons)
    
    answer = inputdlg({'Number of training neurons (or ''all'')'}, 'Nonlinear fit', 1, {'all'});
    if isempty(answer), return; end
    ntrain = str2double(answer{1});
    
end
if strcmpi(answer{1}, 'all')
    
    training_neurons = arrayfun(@(v) setdiff(1:nneurons, v), 1:nneurons, 'uniformoutput', false);
    training_neurons = cat(1, training_neurons{:});
    testing_neurons = [];
    for j = 1:size(training_neurons, 1)
        
        testing_neurons = [testing_neurons; setdiff(1:nneurons, training_neurons(j, :))]; %#ok<AGROW>
        
    end

    
elseif ntrain == 1
    
    training_neurons = [2:nneurons, 1]';  % could randomize order, but for now let's not, for reproducibility
    testing_neurons = [];
    for j = 1:size(training_neurons, 1)
        
        testing_neurons = [testing_neurons; setdiff(1:nneurons, training_neurons(j, :))]; %#ok<AGROW>
        
    end

    
else
    
    assert(ntrain <= nneurons - 1);
    testing_neurons = (1:nneurons)';
    training_neurons = [];
    for j = 1:nneurons
        
        others = setdiff(1:nneurons, j);
        training_neurons = [training_neurons; others(randperm(nneurons - 1, ntrain))];         %#ok<AGROW>
        
    end
    
end


function oedb = assign_params_toseg(oedb, fileinfo, d, n, s, algind, newP, indicatorstring)
fn = fieldnames(newP);
r = fetch_results(oedb, fileinfo, d, n, s, algind);
if ~isfield(newP, 'kd_ex')
    
    newP.kd_ex = [];  %keep empty, don't assign a default value
    
end
if isempty(fieldnames(r.params))
    
    P = sbm.init.assign_default_settings_and_params(struct(), struct(), indicatorstring);
    
else
    
    P = r.params;
    
end

for u = 1:numel(fn)
    
    P.(fn{u}) = newP.(fn{u});
    
end

r.params = orderfields(P);
oedb = oedb_assign_results(oedb, fileinfo, d, n, s, algind, r);