![gpxrdpy](./docs/gpxrdpy_banner_light.svg#gh-light-mode-only)
![gpxrdpy](./docs/gpxrdpy_banner_dark.svg#gh-dark-mode-only)

Python wrapper for PXRD pattern calculation based on pyobjcryst. There are four main usages:
* calculate PXRD pattern based on cif file - `calc`
* compare two PXRD patterns based on .tsv files - `comp`
* fit a background profile - `background`
* average PXRD pattern over the course of a trajectory (.h5) - `average`

Each of these modes can be called using the main script, for which a specific `--help` has been implemented. <br>
Usage:

`gpxrd.py mode [<options>]`
  
Warning: the `average` mode is deprecated in favor of the included gpxrd_average.sh script, as the current backend suffers from some kind of memory error. <br>
Usage:

`gpxrd_average.sh file.h5 nr_snapshots`

The default values for the `run_up_time` and `exp_filename` should be altered within this bash script before execution.

A version of gpxrdpy is also implemented in pyiron at https://github.com/SanderBorgmans/pyiron/tree/hpc_ugent_2020. An example notebook, created as supporting information for [this paper](https://doi.org/10.1039/D3TA00470H), can be found in the notebooks folder. In this notebook, static and dynamic simulations are performed, and provided as input to the gpxrdpy module of pyiron. 

## How to cite 
[1] S. Borgmans, S.M.J. Rogge, J.S. De Vos, C.V. Stevens, P. Van Der Voort, V. Van Speybroeck. Quantifying the Likelihood of Structural Models through a Dynamically Enhanced Powder X-Ray Diffraction Protocol. Angew. Chem. Int. Ed., 60, 16 (2021). https://doi.org/10.1002/anie.202017153

## Requirements
* `pyobjcryst` - Python bindings to ObjCryst++, the Object-Oriented Crystallographic Library (see https://github.com/diffpy/pyobjcryst)
* `numpy` - library for scientific computing in python
* `scipy` - library for scientific computing in python
* `molmod` - collection of molecular modelling tools for python (https://github.com/molmod/molmod)

Additional installation instructions can be found in the [documentation](https://molmod.github.io/gpxrdpy/).
