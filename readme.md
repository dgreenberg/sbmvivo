# sbmvivo
***Disclamer: the API and software in this repository are currently in a pre-alpha phase and subject to change without notice. Backward compatibility might not be maintained.***

This repository contains the sequential binding model (SBM) and tools for model-fitting and AP inference on *in vivo* fluorescence data.

The model and related algorithms are described in our [bioRxiv paper](). Tools for *in vitro* data can be found in [sbmvitro](https://github.com/dgreenberg/sbmvitro).

`sbmvivo` also contains wrappers for several other algorithms and a graphical user interface to compare them over multiple recordings, neurons and datasets.

Most of the code is written in Matlab, with a few performance-critical routines coded in C++. The AP inference routine for the SBM is implemented in both Matlab and CUDA C++ for testing purposes, but computation using the CUDA version on a GPU is essentially mandatory when inferring neuron-specific parameters on all but the smallest datasets.

## Installation
After cloning the repository, open Matlab, change your current directory to the base directory of this repository and run `install_sbmvivo` to update your path. The path changes will be saved if you have write access to your `pathdef.m`.

### Dependencies
Fitting all SBM parameters to a dataset with known AP times requires the optimization toolbox, while AP inference and parameter inference from fluorescence alone do not.

CUDA algorithms require [CUB](https://nvlabs.github.io/cub/). Pre-compiled MEX wrappers for CUDA algorithms should work out of the box on most systems, but otherwise a readme and makefile are available for recompilation. Currently they have only been compiled and tested on Linux, though non-CUDA MEX files work on windows as well.

## Command line tools
SBM-based inference of AP times and neuron-specific parameters from fluorescence alone, when global parameters are known, is available through the command
```matlab
sbm.infer(fluorescence, image_times, indicator)
```
see `help(sbm.infer)` for more details. Currently, the valid choices for the `indicator` string input are 'gcamp6s', 'gcamp6f' and 'ogb1-am'. Default global parameters will be used depending on the indicator, if available. Alternatively, you can supply your own parameters with
```matlab
sbm.infer(fluorescence, image_times, indicator, params)
```
Any parameters specified in `params` will override indicator-specific defaults. To infer AP times only *without inferring neuron-specific parameters*, use
```matlab
opts.parameterestimation = 'none';
sbm.infer(fluorescence, image_times, indicator, params, opts)
```

To fit all global and neuron-specific parameters to fluorescence with known AP times, use the function `sbm.fit.globalfit()` or its wrapper `gui/oedb_nonlinfit`.

Additional API documentation will be added in the near future.

## GUI
After installation, the graphical user interface can started with the command `oedatabrowser`. If the option "keep data in memory" is not selected, a temporary directory will be created to store the data and results, which can be useful when these do not fit in memory.

For SBM-based AP inference *without inferring neuron-specific parameters*, click the "filter/smooth" button.

Datasets can be loaded into `oedatabrowser` automatically on startup by placing them in the data subdirectory of this repository, or can be manually imported.

Pressing **ctrl + D** in `oedatabrowser` saves several variables describing the currently selected data and AP inference results to the base Matlab workspace.

## Data file formats
Datasets should be saved in `.mat` files containing a variable named `oerec`. The data format is described in `structs/empty_oerec` and other files referenced therein. An example can be found in the supplementary data of the [bioRxiv paper]().

Libraries containing multiple datasets, parameter sets, algorithm settings and inferred APs are saved and loaded by `oedatabrowser` and the command line tools as `.odb` files, which are actually zip files containing `.mat` files.

## Adding additional algorithms
Additional methods can be added to the package by putting a `.m` file beginning with `apdet_` in the algorithms folder. This file should contain a function with the same name as the file, which takes three inputs `(oerec, params, opts)` and return two outputs `(results, alginfo)`.

In the near future, it will also be possible to add training routines for parameter fitting and cross-validation for new algorithms.

## GPU debugging
The GPU-based inference method can also be run as a standalone application for easier debugging. The Matlab-based command-line tool `sbm.toraw()` and GUI can export an inference problem to a file that can be read by this application. This process will be further-documented in the near future.
