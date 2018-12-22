# sbmvivo
***Disclamer: the API and software in this repository are currently in a pre-alpha phase and subject to change without notice. Backward compatibility might not be maintained.***

This repository contains the sequential binding model (SBM) and tools for model-fitting and AP inference on *in vivo* fluorescence data.

The model and related algorithms are described in our [bioRxiv paper](https://www.biorxiv.org/content/early/2018/11/29/479055). Tools for *in vitro* data can be found in [sbmvitro](https://github.com/dgreenberg/sbmvitro).

`sbmvivo` also contains wrappers for several other algorithms and a graphical user interface to compare them over multiple recordings, neurons and datasets.

Most of the code is written in Matlab, with a few performance-critical routines coded in C++. The AP inference routine for the SBM is implemented in both Matlab and CUDA C++ for testing purposes, but computation using the CUDA version on a GPU is essentially mandatory when inferring neuron-specific parameters on all but the smallest datasets. A GPU with CUDA compute capability 3.0 or higher is required.

## Installation
***Warning: do not simply add the entire repository with all subfolders to the Matlab path!***

After cloning the repository, open Matlab, change your current directory to the base directory of this repository and run `install_sbmvivo` to update your path. The path changes will be saved if you have write access to your `pathdef.m`.

### Dependencies
Fitting all SBM parameters to a dataset with known AP times requires the optimization toolbox, while AP inference and parameter inference from fluorescence alone do not.

CUDA algorithms require [CUB](https://nvlabs.github.io/cub/) for compilation only.

### Precompiled CUDA MEX files
Precompiled MEX wrappers for CUDA algorithms are provided for the following setup:

* Linux
* CUDA 9.1 release 9.1.85
* Nvidia driver version 390.48
* GPU: Tesla K40, Tesla K20c, GeForce GTX 1080 or GeForce RTX 2080
* Matlab R2018a or R2018b

They may also work with other configurations but have not been widely tested. Changing the graphics card or Matlab version probably won't make a difference, while changing the CUDA version may require recompilation. If the precompiled MEX files do not work out of the box on your system, you will have to recompile them (see below).

## Command line tools
SBM-based inference of AP times and neuron-specific parameters from fluorescence alone, when global parameters are known, is available through the command
```matlab
aptimes = sbm.infer(fluorescence, image_times, indicator);
```
see `help(sbm.infer)` for more details.

As an example of how to prepare the inputs to `sbm.infer`, suppose we want to run the SBM on all data for the third neuron of the dataset provided with the arXiv paper. Then we would use
```matlab
load('in vivo imaging with AP times.mat');
fluorescence = {oerec(3).data.f};
image_times = {oerec(3).data.t};
indicator = 'gcamp6s';
```

Currently, the only valid choices for the `indicator` string input is 'gcamp6s', but this will change in the near future. Default global parameters will be used depending on the indicator, if available. Alternatively, you can supply your own parameters with
```matlab
aptimes = sbm.infer(fluorescence, image_times, indicator, params)
```
Any parameters specified in `params` will override indicator-specific defaults. To infer AP times only *without inferring neuron-specific parameters*, use
```matlab
opts.parameterestimation = 'none';
aptimes = sbm.infer(fluorescence, image_times, indicator, params, opts);
```

To use the CPU instead of the GPU for computation, use:
```matlab
opts.usegpu = false;
```
However, GPU computation is essential for estimation of neuron-specific parameters from fluorescence data alone in a reasonable time, and useful in other cases as well.

To fit all global and neuron-specific parameters to fluorescence with known AP times, use the function `sbm.fit.globalfit()` or its wrapper `gui/oedb_nonlinfit`.

Additional API documentation will be added in the near future.

## GUI
After installation, the graphical user interface can be started with the command `oedatabrowser`. If the option "keep data in memory" is not selected, a temporary directory will be created to store the data and results, which can be useful when these do not fit in memory.

For SBM-based AP inference *without inferring neuron-specific parameters*, click the "filter/smooth" button.

Datasets can be loaded into `oedatabrowser` automatically on startup by placing them in the data subdirectory of this repository, or can be manually imported.

To disable GPU computation, uncheck the menu item "Options | SBM | GPU computation".

Pressing **ctrl + D** in `oedatabrowser` saves several variables describing the currently selected data and AP inference results to the base Matlab workspace.

## Data file formats
Datasets should be saved in `.mat` files containing a variable named `oerec`. The data format is described in `structs/empty_oerec` and other files referenced therein. An example can be found in the supplementary data of the [bioRxiv paper](https://www.biorxiv.org/content/early/2018/11/29/479055).

Libraries containing multiple datasets, parameter sets, algorithm settings and inferred APs are saved and loaded by `oedatabrowser` and the command line tools as `.odb` files, which are actually zip files containing `.mat` files.

## Adding additional algorithms
Additional methods can be added to the package by putting a `.m` file beginning with `apdet_` in the algorithms folder. This file should contain a function with the same name as the file, which takes three inputs `(oerec, params, opts)` and return two outputs `(results, alginfo)`.

In the near future, it will also be possible to add training routines for parameter fitting and cross-validation for new algorithms.

## GPU debugging
The GPU-based inference method can also be run as a standalone application for easier debugging. The Matlab-based command-line tool `sbm.toraw()` and GUI can export an inference problem to a file that can be read by this application. This process will be further-documented in the near future.

## MEX file compilation
Non-CUDA MEX files can be compiled directly using MEX. Some have comments in the file indicating how they should be compiled, for an example see `fspikecountsfit_5state_odesolver_backwardeuler_mainloop.cpp` in the subdirectory `mex/sbm_fit`.

CUDA MEX files should be compiled using a makefile, which is found in the subdirectory `mex/sbm_cuda`. So far, they have only been compiled on linux. If you encounter a compilation or runtime error, you can create an issue on this repository to ask for help. In that case, please make sure to include the following:
* Information about your OS and Matlab version
* Output of the command nvcc --version
* Output of the command nvidia-smi
