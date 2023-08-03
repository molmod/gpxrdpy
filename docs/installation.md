---
hide:
  - toc
  - navigation
---

## Setup

The easiest way to install gpxrdpy, is by first installing the pyobjcryst through conda:

```console
  $ conda create -n gpxrdpy
  $ conda activate gpxrdpy
  $ conda install -c conda-forge pyobjcryst
```

Afterward, the gpxdpy package can be installed through pip:

```console
  $ pip install --upgrade git+https://github.com/SanderBorgmans/gpxrdpy.git
```

## Checking integrity of installation

To quickly check whether the installation went smoothly, you can run the gpxrdpy-test.py script as follows:

```console
  $ gpxrdpy-test.py
```

which computes the static PXRD pattern of COF-5. This will generate a `.fhkl` file with the intensities for each hkl index, and a `.tsv` file with the integrated pattern. 