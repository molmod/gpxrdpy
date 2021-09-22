#!/usr/bin/env python

import sys,h5py,os,glob
from optparse import OptionParser
import numpy as np

from molmod.units import *
from gpxrdpy.gpxrd_utils import *


if __name__ == "__main__":
    if (len(sys.argv) < 3):
        print("Usage: gpxrd.py mode [<options>]")
        print("modes: calc, comp, average, sp, prepare, post, background")
        print("Use the --help option per mode for more info.")
        sys.exit(1)

    if not sys.argv[1] in ["calc", "comp", "average", "sp", "prepare", "post","background"]:
        print("Invalid mode: {}".format(sys.argv[1]))
        print("Usage: gpxrd mode [<options>]")
        print("modes: calc, comp, average, sp, prepare, post, background")
        print("Use the --help option per mode for more info.")
        sys.exit(1)

    # Calc mode
    if sys.argv[1]=="calc":
        parser = OptionParser(usage="Usage: %prog calc <ciffile> [<options>]")
        parser.add_option("--wavelength",
                  action="store", type="float", dest="wavelength", help="set wavelength (Angstroms) [default: %default]", default=1.54056)
        parser.add_option("--peakwidth",
                  action="store", type="float", dest="peakwidth", help="set peak width (degrees) [default: %default]", default=0.14)
        parser.add_option("--numpoints",
                  action="store", type="int", dest="numpoints", help="set number of points [default: %default]", default=1000)
        parser.add_option("--max2theta",
                  action="store", type="float", dest="max2theta", help="set max 2*theta (degrees) [default: %default]", default=50)
        parser.add_option("--obspattern",
                  action="store", type="str", dest="obspattern", help="set file name to read observed pattern [default: %default]", default=None)
        parser.add_option("--plot",
                  action="store", type="str", dest="plot", help="plot the data (to file if specified) [default: %default]", default='False')
        parser.add_option("--check_peaks",
                  action="store_false", dest="check_peaks", help="check whether there are significant peaks outside the obspattern range [default: %default]", default=True)
        parser.add_option("--detail",
                  action="store_true", dest="detail", help="print detailed Fhkl info per hkl [default: %default]", default=False)

        (options, args) = parser.parse_args(sys.argv[2:])

        # Process the arguments
        assert len(args)==1
        filename = args[0]

        # Basic argument checks
        if not filename.split('.')[-1] == 'cif':
            print("The file extension is faulty \n")
            print("Usage: gpxrd calc <cifFile> [<options>]\n")
            print("Use --help for more info.\n")
            sys.exit(1)

        wavelength  = options.wavelength
        peakwidth   = options.peakwidth * deg #this is the FWHM
        numpoints   = options.numpoints
        max2theta   = options.max2theta * deg
        obspattern  = options.obspattern
        plot        = options.plot
        check_peaks = options.check_peaks
        detail      = options.detail

        # Calculate the PXRD pattern and save the output
        crystal = create_crystal(filename)
        calculate(filename, crystal, wavelength, peakwidth, numpoints, max2theta, obspattern, plot, check_peaks, detail)


    # Comp mode
    if sys.argv[1]=="comp":
        parser = OptionParser(usage="Usage: %prog comp <pattern1-file> <pattern2-file> [<options>]\n"
                                    "The pattern files should contain two tab separated columns, 2theta and intensity.\n"
                                    "Internally the first pattern is expected to be the reference pattern.\n")
        parser.add_option("--plot",
                  action="store", type="str", dest="plot", help="plot the data (to file if specified) [default: %default]", default='False')
        parser.add_option("--no_scale",
                  action="store_false", dest="scale", help="set this flag if you don't want to scale the patterns", default=True)
        parser.add_option("--scale_max",
                  action="store_true", dest="scale_max", help="set this flag if you want want a fixed max intensity", default=False)

        (options, args) = parser.parse_args(sys.argv[2:])

        # Process the arguments
        assert len(args)==2
        obspattern, calcpattern = args
        plot  = options.plot
        scale = options.scale
        scale_max = options.scale_max

        # Compare the two patterns and save the output
        compare(obspattern, calcpattern, plot, scale=scale, scale_max=scale_max)

    # average mode
    if sys.argv[1]=="average":
        raise DeprecatedError("Please use a combination of the prepare and sp modes, because memory errors arise!")
        parser = OptionParser(usage="Usage: %prog average <traj.h5> <pattern-file> [<options>]\n"
                                    "An average PXRD pattern will be generated from snapshots taken from the trajectory file (in .h5 format).\n")
        parser.add_option("--runuptime",
                  action="store", type="int", dest="runuptime", help="define run up time in number of steps [default: %default]", default=1000)
        parser.add_option("--snapshots",
                  action="store", type="int", dest="snapshots", help="set number of snapshots [default: %default]", default=1000)
        parser.add_option("--wavelength",
                  action="store", type="float", dest="wavelength", help="set wavelength (Angstroms) [default: %default]", default=1.54056)
        parser.add_option("--peakwidth",
                  action="store", type="float", dest="peakwidth", help="set peak width (degrees) [default: %default]", default=0.14)
        parser.add_option("--plot",
                  action="store", type="str", dest="plot", help="plot the data (to a file if you specify a filename) [default: %default]", default='False')
        parser.add_option("--check_peaks",
                  action="store_false", dest="check_peaks", help="check whether there are significant peaks outside the obspattern range [default: %default]", default=True)

        (options, args) = parser.parse_args(sys.argv[2:])

        # Process the arguments
        assert len(args)==2
        filename   = args[0]
        obspattern = args[1]

        # Basic argument checks
        if not filename.split('.')[-1] == 'h5':
            print("The file extension is faulty \n")
            print("Usage: gpxrd average <traj.h5> <pattern-file> [<options>]\n")
            print("Use --help for more info.\n")
            sys.exit(1)

        runuptime   = options.runuptime
        snapshots   = options.snapshots
        wavelength  = options.wavelength
        peakwidth   = options.peakwidth * deg
        plot        = options.plot
        check_peaks = options.check_peaks

        patterns = []
        try:
            os.mkdir('cifs')
        except FileExistsError:
            print('Attention, the "cifs" directory already exists!')

        with h5py.File(filename, 'r') as traj:
            idx = snapshots2cifs(traj,runuptime,snapshots)

        # Create crystals
        crystals = create_crystals(idx)

        # Calculate patterns
        patterns = np.zeros(0)
        warnings = []
        for n,crystal in enumerate(crystals):
            print(n)
            ttheta, iobs, pattern, warning  = calculate_pattern(crystal, wavelength, peakwidth, obspattern, check_peaks)
            warnings.append(warning)
            if len(patterns)==0:
                patterns = np.concatenate((patterns,pattern))
            else:
                patterns = np.vstack((patterns,pattern))

        #ttheta, iobs, patterns = calculate_patterns(crystals, wavelength, peakwidth, obspattern, check_peaks)
        icalc = np.average(patterns,axis=0)
        scalefactor = FitScaleFactorForRw(iobs,icalc,iobs.max()/icalc.max())

        statistical_comparison(ttheta,iobs,icalc*scalefactor)
        plot_data(ttheta,iobs,icalc*scalefactor,plot)



    # prepare mode (executed before any sp calculation)
    # This mode should be deprecated soon if the code can handle more than 1 pattern at the same time in memory
    if sys.argv[1]=="prepare":
        parser = OptionParser(usage="Usage: %prog prepare <traj.h5> [<options>]\n"
                                    "Prepares the snapshots for average PXRD pattern taken from the trajectory file (in .h5 format).\n")
        parser.add_option("--runuptime",
                  action="store", type="int", dest="runuptime", help="define run up time in number of steps [default: %default]", default=1000)
        parser.add_option("--snapshots",
                  action="store", type="int", dest="snapshots", help="set number of snapshots [default: %default]", default=1000)

        (options, args) = parser.parse_args(sys.argv[2:])

        # Process the arguments
        assert len(args)==1
        filename   = args[0]

        # Basic argument checks
        if not filename.split('.')[-1] == 'h5':
            print("The file extension is faulty \n")
            print("Usage: gpxrd prepare <traj.h5> [<options>]\n")
            print("Use --help for more info.\n")
            sys.exit(1)

        runuptime  = options.runuptime
        snapshots  = options.snapshots

        try:
            os.mkdir('cifs')
        except FileExistsError:
            print('Attention, the "cifs" directory already exists!')


        traj = h5py.File(filename, 'r')
        idx = snapshots2cifs(traj,runuptime,snapshots)
        traj.close()

    # sp mode
    # This mode can be deprecated if the underlying C code can handle more than 1 pattern at the same time in memory
    if sys.argv[1]=="sp":
        parser = OptionParser(usage="Usage: %prog sp <idx> <pattern-file> [<options>]\n"
                                    "A PXRD pattern will be generated from this snapshots and stored.\n")
        parser.add_option("--wavelength",
                  action="store", type="float", dest="wavelength", help="set wavelength (Angstroms) [default: %default]", default=1.54056)
        parser.add_option("--peakwidth",
                  action="store", type="float", dest="peakwidth", help="set peak width (degrees) [default: %default]", default=0.14)
        parser.add_option("--check_peaks",
                  action="store_false", dest="check_peaks", help="check whether there are significant peaks outside the obspattern range [default: %default]", default=True)

        (options, args) = parser.parse_args(sys.argv[2:])

        # Process the arguments
        assert len(args)==2
        idx        = int(args[0])
        obspattern = args[1]

        wavelength  = options.wavelength
        peakwidth   = options.peakwidth * deg
        check_peaks = options.check_peaks

        try:
            os.mkdir('frames')
        except FileExistsError:
            print('Attention, the "frames" directory already exists!')

        # Create crystals
        crystal = create_crystals([idx])[0]

        # Calculate patterns
        ttheta, iobs, icalc, warning, fttheta, ficalc = calculate_pattern(crystal, wavelength, peakwidth, obspattern, check_peaks)
        # Do not scale the patterns, only scale the fully averaged one
        #ficalc = ficalc * iobs.max()/icalc.max()
        #icalc = icalc * iobs.max()/icalc.max()
        sp_output(ttheta,icalc,idx,warning,fttheta,ficalc)

    # post mode
    # This mode should be deprecated soon if the code can handle more than 1 pattern at the same time in memory
    if sys.argv[1]=="post":
        parser = OptionParser(usage="Usage: %prog post <pattern-file> [<options>]\n"
                                    "An average PXRD pattern will be generated from the stored snapshots and compared.\n")
        parser.add_option("--plot",
                  action="store", type="str", dest="plot", help="plot the data (to a file if you specify a filename) [default: %default]", default='False')

        (options, args) = parser.parse_args(sys.argv[2:])

        # Process the arguments
        assert len(args)==1
        obspattern = args[0]

        plot = options.plot

        # Average over all stored patterns and store it
        warnings = []
        pattern = np.load('frames/norm_0.npy')
        warning,pattern = filter_nan(pattern)
        warnings.append(warning)
        nframes = len(glob.glob('frames/norm_*.npy'))
        for i in range(1,nframes):
            data = np.load('frames/norm_{}.npy'.format(i))
            warning,data = filter_nan(data)
            warnings.append(warning)
            pattern[:,1] += data[:,1]
        pattern[:,1] /= nframes


        fpattern = np.load('frames/full_0.npy')
        nframes = len(glob.glob('frames/full_*.npy'))
        for i in range(1,nframes):
            fpattern[:,1] += np.load('frames/full_{}.npy'.format(i))[:,1]
        fpattern[:,1] /= nframes


        with open('avg.dat','w') as f:
            for i in range(pattern.shape[0]):
                f.write("{:4.3f}\t{:10.8f}\n".format(pattern[i,0],pattern[i,1]))

        with open('favg.dat','w') as f:
            for i in range(fpattern.shape[0]):
                f.write("{:4.3f}\t{:10.8f}\n".format(fpattern[i,0],fpattern[i,1]))

        # Calculate patterns
        compare('avg.dat', obspattern, plot, warning=any(warnings), full_data=fpattern)

    if sys.argv[1]=="background":
        parser = OptionParser(usage="Usage: %prog background <pattern-file> [<options>]\n"
                                    "The background will be estimated.\n")
        parser.add_option("--locs",
                  action="callback", type="string", dest="locs", callback=callback, help="node locations in degrees, comma separated list")
        parser.add_option("--range",
                  action="callback", type="string", dest="range", callback=callback, help="range for background points, comma separated list")
        parser.add_option("--bkg_points",
                  action="store", type="int", dest="bkg_points", help="set number of background points in fit [default: %default]", default=10)
        parser.add_option("--uniform",
                  action="store_true", dest="uniform", help="no longer adapt the uniform grid for interpolation to all minima [default: %default]", default=False)
        parser.add_option("--wavelength",
                  action="store", type="float", dest="wavelength", help="set wavelength (Angstroms) [default: %default]", default=1.54056)
        parser.add_option("--peakwidth",
                  action="store", type="float", dest="peakwidth", help="set peak width (degrees) [default: %default]", default=0.14)
        parser.add_option("--max2theta",
                  action="store", type="float", dest="max2theta", help="set max 2*theta (degrees) [default: %default]", default=50)
        parser.add_option("--plot",
                  action="store", type="str", dest="plot", help="plot the data (to a file if you specify a filename) [default: %default]", default='False')

        (options, args) = parser.parse_args(sys.argv[2:])

        # Process the arguments
        assert len(args)==1
        obspattern = args[0]

        range      = options.range
        locs       = options.locs
        bkg_points = options.bkg_points
        uniform    = options.uniform
        wavelength = options.wavelength
        peakwidth  = options.peakwidth * deg
        max2theta  = options.max2theta * deg

        plot = options.plot
        background(obspattern, peakwidth, bkg_points, locs, range, uniform, plot)
