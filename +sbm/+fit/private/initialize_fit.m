function [simulation_function, p, parameter_names, statenames, fixed_params, params_to_print] = initialize_fit(opts)
%initialize parameters and get additional info about the model
statenames = model_statenames(opts);

p = struct( ...
    'R', 0.05 ...  % ratio of background to neuronal fluorescence, reduces DF/F0
    ,'spiketshift', 0 ...
    ,'concentration_scale_to_nM', 1000 ... we use micromolar units
    ,'c0', 0.05 ...
    ,'fratio', 45 ... fluorescence ratio of calcium saturated to calcium free indicator, minus 1
    ,'S', 20 ... indicator concentration
    );

fixed_params = {'spiketshift' 'concentration_scale_to_nM' 'cafree_perAP_range_uM'};

if opts.adjust_rate_constants
    
    fixed_params = [fixed_params {'c0'}];
    
else
    
    fixed_params = [fixed_params {'koff' 'kon' 'k50' 'tau_mm' 'k50_first' 'dk50'}];
    
end

p.fratio = 45; % fluorescence ratio of calcium saturated to calcium free indicator, minus 1

if opts.nbindingsteps > 1
    
    p.brightness_ratio = 0.5;
    if opts.nbindingsteps > 2
        
        p.brightness_ratio = [0.01 p.brightness_ratio];
        
        if opts.nbindingsteps > 3
            
            p.brightness_ratio = [0.9 * ones(1, opts.nbindingsteps - 3) p.brightness_ratio];
            
        end
        
    end
    
end

%get the size of the AP so we can normalize other params by it
if opts.nbuffers == 0
    
    p.A = 0.5; %uM
    
else
    
    p.A = 5; %uM
    
end

if opts.use_kd_ex
    
    p.kd_ex = ones(1, opts.n_ex) * 0.8; %uM
    
end

%here's a formula we DON'T use to initialize the extrusion time constant:
%
%cbr_ogb = 20 * 0.2 / (0.2 + 0.05) ^ 2;  % 64 = calcium binding ratio of OGB at 20 uM from Stosiek et al. 2003, assuming .05 uM c0 and 20 uM OGB1 Kd
%gamma = (1 + target_buffering_capacity + cbr_ogb) / 0.3; % formula from Helmchen et al. 1996. empirical decay constant of 300 ms from Stosiek et al. 2003
%
%Why not?
%Because the equilibrium approximation is not quite true for either
%endogenous buffers or indicators, the apparent gamma will be higher
%than the true one as plenty of calcium ions will be extruded before
%they can be bound.

if opts.nbuffers == 0
    
    p.tau_ex = 0.3; %calcium extrusion time constant will seem a lot slower when we don't model buffers
    
elseif opts.use_kd_ex && opts.n_ex > 1
    
    p.tau_ex = linspace(0.005, 0.25, opts.n_ex);
    
else
    
    p.tau_ex = 1 / 50;
    
end

%buffering:
switch opts.nbuffers
    case 0
        
        [p.kd_B, p.tau_B, p.Btot] = deal([]);
        fixed_params = [fixed_params {'kd_B' 'tau_B' 'Btot'}];
        
    case 1
        
        %single site approximation of Calbindin, Naegerl et al. 2000
        p.kd_B = 0.158; %uM
        p.tau_B = 0.2671;
        p.Btot = 1;
        
    case 2
        
        p.kd_B = [0.5 0.5]; %uM
        p.tau_B = [15e-3 1];
        p.Btot = [20 20];
        
        %calbindin 3:1 model, Naegerl et al. 2000
        %p.kd_B = [0.513 0.175]; %uM
        %p.tau_B = [0.0253 0.4396];
        %p.Btot = [1 3]; %3:1 binding site ratio
    otherwise
        
        p.kd_B = 0.5 * ones(1, opts.nbuffers);  %uM
        p.tau_B = linspace(15e-3, 1.5, opts.nbuffers);
        p.Btot = ones(1, opts.nbuffers);
        
end

%set buffer concentrations to achieve the desired buffering capacity at baseline calcium
if opts.nbuffers > 0
    
    fast_buffers = p.tau_B < opts.fast_buffer_threshold;
    if ~any(fast_buffers)
        
        warning('no fast buffer available');
        fast_buffers = p.tau_B == min(p.tau_B);
        
    end
    p.Btot = p.Btot * opts.target_buffering_capacity / sum(p.Btot(fast_buffers) .* p.kd_B(fast_buffers) ./ (p.kd_B(fast_buffers) + p.c0).^2);
    
end

switch opts.modelname
    case {'5sb'}
        
        if ~opts.invitroinit
            
            %             mink50_init = 0.25;
            %             maxk50_init = 0.425;
            %             k50 = linspace(mink50_init, maxk50_init, opts.nbindingsteps);
            k50 = [.325 .37  .875 .88]; %gcamp6s 8sub
            tau_init = 0.1;
            
            %k50 = [0.7500 0.7706    2.2706    3.7706]; % gcamp6f 8sub
            %tau_init = 0.04;
            
            tau = ones(1, opts.nbindingsteps) * tau_init;
            
            if opts.usek50
                
                p.k50_first = k50(1);
                p.dk50 = diff(k50);
                p.tau_mm = tau;
                
            else
                
                [koff, kon] = k50tau2koffkon(k50, tau);
                p.koff = koff;
                p.kon =  kon;
                
            end
            
        else
            
            if true
                
                %from fit50par 200518.pickle
                k50 = cumsum([0.26686,  0.054469,  0.011218,  0.070457]);
                tau_mm = [0.35669 0.0021059 0.080566 0.063355];                                 
                
            elseif false
                
                %from fit50par_040118_SWdata.pickle
                k50 = [0.25506565,  0.25756565,  0.34268284,  0.39195712];
                tau_mm = [0.00227183,  0.00996189,  0.34782586,  0.04774549];                                 
                
            else
                
                %sun paper's C-terminus before N-terminus concept with faas et al. rate constants
                koff = [2.6e3 6.5 1.6e5 2.2e4]; % s^-1
                kon = [84 25 770 3.2e4]; %s^-1 (uM)^-1
                [k50, tau_mm] = koffkon2k50tau(koff, kon);
                
            end
            
            if opts.usek50
                
                p.k50_first = k50(1);
                p.dk50 = diff(k50);
                p.tau_mm = tau_mm;
                
            else
                
                [p.koff, p.kon] = k50tau2koffkon(k50, tau_mm);
                
            end
            
        end
        
        if opts.use_kd_ex
            
            simulation_function = @fspikecountsfit_5stateb_odesolver_multi_kdex;
            
        else
            
            simulation_function = @fspikecountsfit_5stateb_odesolver_multi;
            
        end
        
    otherwise
        error('unrecognized model: %s', modelname);
end
p = orderfields(p);

parameter_names = fieldnames(p);

params_to_print = parameter_names;
if ismember('fratio', params_to_print)
    
    params_to_print = [setdiff(params_to_print(:), {'fratio' 'brightness_ratio'}'); {'dbrightness'}];
    
end
if opts.usek50
    
    params_to_print = [setdiff(params_to_print(:), {'k50_first' 'dk50'}'); {'k50' 'koff' 'kon'}'];  % print both representations
    
    
end
if opts.autoR
    
    parameter_names = setdiff(parameter_names, {'R'}); %we're doing it automatically
    params_to_print = setdiff(params_to_print, {'R'}); %will be printed anyway, but not as a parameter
    
end


function statenames = model_statenames(opts)
statenames = {'Ca_{free}' 'G'};
for k = 1:opts.nbindingsteps
    
    if k == 1
        
        statenames{1, end + 1} = 'CaG';
        
    else
        
        statenames{1, end + 1} = sprintf('Ca_%dG', k);
        
    end
    
end
for k = 1:opts.nbuffers
    
    statenames{1, end + 1} = sprintf('B_%d', k);
    statenames{1, end + 1} = sprintf('CaB_%d', k);
    
end