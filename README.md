# gpxrdpy
Python wrapper for PXRD pattern calculation based on pyobjcryst. There are four main usages:
* calculate PXRD pattern based on cif file - `calc`
* compare two PXRD patterns based on .tsv files - `comp`
* fit a background profile - `background`
* average PXRD pattern over the course of a trajectory (.h5) - `average`

Each of these modes can be called using the main script, for which a specific `--help` has been implemented.
Usage:

`gpxrd.py mode [<options>]`
  
Warning: the `average` mode is deprecated in favor of the included gpxrd_average.sh script, as the current backend suffers from some kind of memory error.
Usage:

`gpxrd_average.sh file.h5 no_snapshots`

The default values for the `run_up_time` and `exp_filename` should be altered within this bash script before execution.

## Requirements
* `pyobjcryst` - Python bindings to ObjCryst++, the Object-Oriented Crystallographic Library (see https://github.com/diffpy/pyobjcryst)
* `numpy` - library for scientific computing in python
* `scipy` - library for scientific computing in python
* `molmod` - collection of molecular modelling tools for python (https://github.com/molmod/molmod)
