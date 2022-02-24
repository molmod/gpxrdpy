import os
import pyobjcryst
import numpy as np
import matplotlib.pyplot as plt
from pyobjcryst.crystal import *
from pyobjcryst.scatteringpower import ScatteringPowerAtom
from pyobjcryst.atom import Atom
from pyobjcryst.polyhedron import MakeTetrahedron
from pyobjcryst.powderpattern import *
from pyobjcryst.radiation import RadiationType
from pyobjcryst.indexing import *
from pyobjcryst.molecule import *
from pyobjcryst.globaloptim import MonteCarlo
from pyobjcryst.io import xml_cryst_file_save_global

from molmod.units import *

px = PowderPattern()
px.SetWavelength("Cu")
px.SetPowderPatternX(np.linspace(0, 50*deg, 1000))
px.SetPowderPatternObs(np.ones(1000)) # Use fake unit intensity observed pattern since no observed pattern

c = CreateCrystalFromCIF('0.cif')

diffData = px.AddPowderPatternDiffraction(c)

px.plot(hkl=True)
