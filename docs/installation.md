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

Afterward, the gpxdpy package can be installed through pip, which will install all other dependencies:

```console
  $ pip install --upgrade git+https://github.com/SanderBorgmans/gpxrdpy.git
```

## Checking integrity of installation

To quickly check whether the installation went smoothly, you can run the gpxrdpy-test.py script as follows:

```console
  $ gpxrd_test.py
```

which computes the static PXRD pattern of COF-5 and compares it to a reference pattern. This will generate a `.fhkl` file with the intensities for each hkl index, and two `.dat` files with the integrated diffraction intensities as a function of the Bragg angle (2θ), and the scattering vector (q). Finally, the comparison between the calculated and the reference pattern is saved to a `xrd.pdf` file, whereas the statistical comparison is printed to the terminal:

```console
  Imported powder pattern: 2351 points, 2theta=  3.000 ->  50.000, step= 0.020
  Guess = 3.725577310703468e-09, fit = 7.3208687727899055e-09

  Statistical comparison
  ------------------------------------------------------------------
  Quantity 		 |     Full    | LA (min-10) | HA (10-max)
  ------------------------------------------------------------------
  R factor (abs) 		 | 0.942803315 | 0.845442566 | 0.985120011
  R factor (squared)	 | 0.924069988 | 0.839995526 | 0.986114068
  Weighted R factor 	 | 0.953021538 | 0.873615415 | 0.985541226
  Similarity index 	 | 0.388794850 | 0.558217914 | 0.476000931

```
