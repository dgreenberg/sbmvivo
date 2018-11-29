function opts = defaultopts(modelname)
%fixme streameline param bound info system

opts.verbose = 10;

%model structure
opts.nbindingsteps = 4;
opts.nbuffers = 2;
opts.use_kd_ex = false;
opts.n_ex = 1;
opts.neuronspecific = {'R' 'S'}; %different for each neuron. generally this is R, S
opts.sessionspecific = {}; %can change for each session, which is a collection of segments (generally a file)
opts.focalplanespecific = {}; %changes when the focal plane shifts
%tell the solver which parameters don't change the states
opts.obsvars = {'R' 'dbrightness' 'brightness_ratio' 'fratio'};

%rate constant info
opts.adjust_rate_constants = true;
opts.max_rate_constant_adjustment = inf;

%initialization info
opts.invitroinit = false;
opts.target_buffering_capacity = 125;
opts.fast_buffer_threshold = 0.3;
opts.mintaufrac = 0.1;
opts.mindiff_init = 1e-5;

%optimization options
opts.autoR = true;
opts.gradreldiff = 5e-7;
opts.gradmindiff = 5e-7;
opts.perispikewin = [0 20];
opts.tau_s0 = 0.5;
opts.tolfun = 1e-9;
opts.tolx = 1e-8;
opts.symgrad = true;
opts.uselog = false;

%(re)parameterization options
opts.usek50 = true;

%parameter bound info. concentrations in micromolar
opts.KDvals = {'kd_B' 'kd_ex'};
opts.minKDval = 0.100;
opts.maxKDval = 100;

if true  %default initialization
    
    opts.min_dk50 = [1e-3 1e-3 1e-3]; %uM
    opts.max_dk50 = 1.5;
    opts.mink50_first = 0.01; %uM
    opts.maxk50_first = 0.75; %uM
    opts.min_tau_mm = 0.001;
    opts.max_tau_mm = 5;
    
else  %restricted sequence etc.
    
    opts.min_dk50 = [0.75    1e-3    1e-3];
    opts.max_dk50 = [1.25    0.06    0.01];
    opts.mink50_first = 0.02; %uM
    opts.maxk50_first = 0.075; %uM
    opts.min_tau_mm = [4.9 0.001 0.05 0.001];
    opts.max_tau_mm = [5 10e-3 150e-3 5e-3];
    opts.min_tau_B = [1e-3 5];
    opts.max_tau_B = [10e-3 20];
    
end

opts.minkon = 0.1; % at 1 uM ligand no forward binding step should have a time constant longer than one second
opts.maxkon = 1e4; % sec ^-1 uM ^-1

opts.minRval = 1e-4;
opts.maxRval = 20;

opts.minc0val = 0.030;
opts.maxc0val = 0.080;

opts.minSval = 0.50;
opts.maxSval = 300;

opts.mintauval = 0.00025;
opts.maxtauval = 30;

opts.max_tau_ex = 0.3;
opts.min_tau_ex = 2.5e-3;
    
opts.minAval = 0.01;
opts.maxAval = 100;
    
opts.minfratio = 10;
opts.maxfratio = 80;