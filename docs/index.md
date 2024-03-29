---
hide:
  - toc
  - navigation
---

![gpxrdpy](./gpxrdpy_banner_light.svg)

Python wrapper for PXRD pattern calculation based on [pyobjcryst](https://github.com/diffpy/pyobjcryst). There are four main usages:

* calculate PXRD pattern based on cif file - `calc`
* compare two PXRD patterns based on .tsv files - `comp`
* fit a background profile - `background`
* average PXRD pattern over the course of a trajectory (.h5) - `average`


## Usage
Each of these modes can be called using the main script, for which a specific `--help` has been implemented. <br>
Usage:

`gpxrd.py mode [<options>]`
  
Warning: the `average` mode is deprecated in favor of the included gpxrd_average.sh script, as the current backend suffers from some kind of memory error. <br>

Usage:

`gpxrd_average.sh file.h5 nr_snapshots`

The default values for the `run_up_time` and `exp_filename` should be altered within this bash script before execution.

A version of gpxrdpy is also implemented in pyiron found in the [ugent pyiron branch](https://github.com/SanderBorgmans/pyiron/tree/hpc_ugent_2020).

## Requirements
* `pyobjcryst` - Python bindings to ObjCryst++, the Object-Oriented Crystallographic Library <br> (see [https://github.com/diffpy/pyobjcryst]())
* `numpy` - library for computing in python
* `scipy` - library for scientific computing in python
* `matplotlib` - library for plotting in python
* `h5py` - library for storing data in python
