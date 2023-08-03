#!/usr/bin/env python
import os

deg = 0.017453292519943295 # conversion of rad to deg


def check_import():
    try:
        import gpxrdpy
    except ModuleNotFoundError:
        print('Could not import gpxrdpy module')

def check_calculation():
    import gpxrdpy
    from gpxrdpy.gpxrd_utils import create_crystal, calculate

    data_path = os.path.join(gpxrdpy.__path__[0], 'data')
    filename = os.path.join(data_path,'COF-5.cif')
    obspattern = os.path.join(data_path,'exp.tsv')

    # Calculate the PXRD pattern and save the output
    # Also provide a reference pattern to perform a comparison analysis

    wavelength  = 1.54056 # angstrom, cu k alpha
    peakwidth   = 0.14 * deg # this is the FWHM
    numpoints   = 1000
    max2theta   = 50 * deg
    obspattern  = obspattern
    plot        = 'False'
    check_peaks = False
    detail      = False # do not print detailed info for each peak
    neutron     = False

    try:
        crystal = create_crystal(filename)
    except Exception as e:
        print('Crystal creation failed')
        print(e)
        return
    
    try:
        calculate(filename, crystal, wavelength, peakwidth, numpoints, max2theta, obspattern, plot, check_peaks, detail, neutron)
    except Exception as e:
        print('PXRD calculation failed')
        print(e)
        return

if __name__ == "__main__":
    print("Testing your gpxrdpy installation...")

    check_import()
    check_calculation()
    