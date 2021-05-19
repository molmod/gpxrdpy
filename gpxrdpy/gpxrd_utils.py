#!/usr/bin/env python
# coding=utf-8

import sys,os,copy,warnings
import numpy as np
import pyobjcryst
import matplotlib.pyplot as pt
from matplotlib.ticker import MaxNLocator
from scipy.optimize import least_squares

from molmod.units import *


##########################################
# Code to redirect stdout
from contextlib import contextmanager

def fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd

@contextmanager
def stdout_redirected(to=os.devnull, stdout=None):
    if stdout is None:
       stdout = sys.stdout

    stdout_fd = fileno(stdout)
    # copy stdout_fd before it is overwritten
    #NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stdout_fd), 'wb') as copied:
        stdout.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(fileno(to), stdout_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(to, 'wb') as to_file:
                os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
        try:
            yield stdout # allow code to be run with the redirected stdout
        finally:
            # restore stdout to its previous value
            #NOTE: dup2 makes stdout_fd inheritable unconditionally
            stdout.flush()
            os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copied

##########################################
# Utility functions for main script
def filter_nan(arr):
    if np.isnan(arr[-1,-1]):
        return True,arr[:-1]
    return False,arr

# taken from https://stackoverflow.com/questions/392041/python-optparse-list
def callback(option, opt, value, parser):
    """
        This function makes a list from a comma separated string
    """
    setattr(parser.values, option.dest, [float(v) for v in value.split(',')])


##########################################
# IO functions

def snapshots2cifs(traj,runuptime,snapshots):
    from molmod.periodic import periodic as ptable # only load this package if necessary

    numbers = traj['/system/numbers'][:]
    coords = traj['/trajectory/pos'][runuptime:]
    rvecs = traj['/trajectory/cell'][runuptime:]

    def get_angle(r0,r1):
        return np.arccos(np.clip((r0*r1).sum(axis=0)/np.linalg.norm(r0)/np.linalg.norm(r1), -1,1))/deg

    idx = np.linspace(0,len(coords)-1,snapshots,dtype=int)
    for n,i in enumerate(idx):
        pos = coords[i]
        rvec = rvecs[i]

        frac = np.dot(pos, np.linalg.inv(rvec))
        alpha = get_angle(rvec[1],rvec[2])
        beta = get_angle(rvec[2],rvec[0])
        gamma = get_angle(rvec[0],rvec[1])

        fn = "cifs/{}.cif".format(n)
        with open(fn,'w') as f:
            f.write("data_h5_{}\n\n".format(i))
            f.write("_cell_length_a\t\t{:.6f}\n".format(np.linalg.norm(rvec[0])/angstrom))
            f.write("_cell_length_b\t\t{:.6f}\n".format(np.linalg.norm(rvec[1])/angstrom))
            f.write("_cell_length_c\t\t{:.6f}\n".format(np.linalg.norm(rvec[2])/angstrom))
            f.write("_cell_angle_alpha\t{:.6f}\n".format(alpha))
            f.write("_cell_angle_beta\t{:.6f}\n".format(beta))
            f.write("_cell_angle_gamma\t{:.6f}\n".format(gamma))
            f.write("_symmetry_space_group_name_H-M\t'P 1'\n_symmetry_Int_Tables_number\t\t1\n\n")
            f.write("loop_\n\t_symmetry_equiv_pos_as_xyz\n\tx,y,z\n\n")
            f.write("loop_\n\t_atom_site_label\n\t_atom_site_type_symbol\n\t_atom_site_fract_x\n\t_atom_site_fract_y\n\t_atom_site_fract_z\n")
            for i in range(len(numbers)):
                s = ptable[numbers[i]].symbol
                f.write("\t{}\t{}\t{:.6f}\t{:.6f}\t{:.6f}\n".format(s+str(i), s, frac[i,0], frac[i,1], frac[i,2]))
    return np.arange(snapshots,dtype=int)

def create_crystal(filename):
    pCryst = pyobjcryst.loadCrystal(filename)
    pCryst.SetUseDynPopCorr(True)
    return pCryst

def create_crystals(idx):
    crystals = []
    for i in idx:
        crystals.append(create_crystal("cifs/{}.cif".format(i)))
    return crystals


##########################################
# Util functions

def R_Factor(icalc,iobs,weighted=False,abs=False):
    # R factor of 0 means a perfect fit
    if weighted:
        wi = [1./iobsi if iobsi>0 else 0. for iobsi in iobs]
        return np.sqrt((wi*(iobs-icalc)*(iobs-icalc)).sum() /(wi*(iobs*iobs)).sum())
    else:
        if abs:
            return (np.abs(iobs-icalc)).sum() / iobs.sum()
        else:
            return np.sqrt(((iobs-icalc)*(iobs-icalc)).sum() /((iobs*iobs)).sum())


def FitScaleFactorForRw(p1,p2,guess):
    """
        Fit scale factor, while keeping observed experimental pattern p2 fixed as reference
    """
    p1 -= np.min(p1)
    p2 -= np.min(p2)
    if np.sum(np.abs(p2)) == 0: return guess # this happens for virtual observed pattern

    def error(x,p1,p2):
        #return R_Factor(p1*x[0],p2,abs=True)
        return R_Factor(p1*x[0],p2,weighted=True)
        #return similarity_index(p1*x[0],p2)

    res = least_squares(error,1.,args=(p1*guess,p2),bounds=(1e-7,np.inf))
    print("Guess = {}, fit = {}".format(guess,guess*res.x[0]))
    return guess*res.x[0]

def similarity_index(p1,p2):
    return (p1*p2).sum()/(np.sqrt((p1*p1).sum()) * np.sqrt((p2*p2).sum()))

def area_between_curves(p1,p2):
    """
        Computes the area between two curves, assuming a uniform grid
    """
    area = 0.
    for i in range(1,len(p1)):
        # calculate area of two triangles and one rectangle
        rect = np.abs(np.min(p1[i-1:i+1]) - np.min(p2[i-1:i+1]))
        triangle1 = np.abs(p1[i] - p1[i-1])/2.
        triangle2 = np.abs(p2[i] - p2[i-1])/2.
        area += rect+triangle1+triangle2
    return area

def statistical_comparison(x,p1,p2,warning=False):
    """
        Statistical comparison, using p2 as reference data
    """
    p1 -= np.min(p1)
    p2 -= np.min(p2)

    if np.sum(np.abs(p2)) == 0: return # this happens for virtual observed pattern

    low_p1 = p1[x<=10*deg]
    low_p2 = p2[x<=10*deg]

    high_p1 = p1[x>=10*deg]
    high_p2 = p2[x>=10*deg]

    print("")
    print("Statistical comparison")
    if warning:
        print("WARNING: the calculated pattern has significant peaks outside the shown range!")
    print("------------------------------------------------------------------")
    print("Quantity \t\t |     Full    | LA (min-10) | HA (10-max)")
    print("------------------------------------------------------------------")
    print("R factor (abs) \t\t | {:10.9f} | {:10.9f} | {:10.9f}".format(R_Factor(p1,p2,abs=True),R_Factor(low_p1,low_p2,abs=True),R_Factor(high_p1,high_p2,abs=True)))
    print("R factor (squared)\t | {:10.9f} | {:10.9f} | {:10.9f}".format(R_Factor(p1,p2),R_Factor(low_p1,low_p2),R_Factor(high_p1,high_p2)))
    print("Weighted R factor \t | {:10.9f} | {:10.9f} | {:10.9f}".format(R_Factor(p1,p2,weighted=True),R_Factor(low_p1,low_p2,weighted=True),R_Factor(high_p1,high_p2,weighted=True)))
    print("Similarity index \t | {:10.9f} | {:10.9f} | {:10.9f}".format(similarity_index(p1,p2),similarity_index(low_p1,low_p2),similarity_index(high_p1,high_p2)))


def plot_data(ttheta,p1,p2,label1='calc',label2='observed',vlines=None,delta=False,warning=False,full_data=None):

    fig = pt.figure()
    ax1 = pt.gca()

    if p2 is None:
        # No reference pattern
        ax1.plot(ttheta/deg,p1,lw=1)
    else:
        # If we are doing a background analysis do not remove anything from p1
        if vlines is None:
            p1 -= np.min(p1)
        p2 -= np.min(p2)

        # Consider the difference for the delta plot
        height = p2.max()-p2.min()
        diff = p2-p1

        # Plot both patterns
        ax1.plot(ttheta/deg,p2,lw=1,label=label2)
        ax1.plot(ttheta/deg,p1,lw=1,label=label1)
        if delta:
            ax1.plot(ttheta/deg,diff-height*0.1,lw=1,label=r'$\Delta$')


    ax1.set_xlabel('2θ (°)')
    ax1.set_ylabel('Intensity (a.u.)')

    if warning:
        ax1.set_title('WARNING: the calculated pattern has significant peaks outside the shown range!')
        if all([data is not None for data in full_data]):
            ax1.plot(full_data[:,0]/deg, full_data[:,1],lw=1,color='gray',alpha=0.5)

    ax1.tick_params(
        axis='y',          # changes apply to the y-axis
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the left edge are off
        right=False,         # ticks along the right edge are off
        labelleft=False) # labels along the left edge are off

    ax1.legend(bbox_to_anchor=(1.1,.5), loc='center left',frameon=False)
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

    lims = pt.xlim()
    ax1.hlines(0,lims[0],lims[1],lw=0.1)
    #ax1.hlines(-height*0.1,lims[0],lims[1],lw=0.1)
    ax1.set_xlim(lims)

    lims = pt.ylim()
    if vlines is not None:
        ax1.vlines(vlines/deg,lims[0],lims[1],lw=0.1)

    if plot.lower() in ['true','yes','y']:
        fig.tight_layout()
        pt.show()
    else:
        path = os.path.splitext(plot)[0] # remove extension
        pt.savefig(path+'.pdf',bbox_inches='tight')
    pt.close()

##########################################
# Actual functions that are executed by main code

def calculate_pattern(crystal, wavelength, peakwidth, pattern, check_peaks):
    # Basic argument check
    if (wavelength < 1.e-6):
        print("Warning: Wavelength is too small. A crash may be iminent.\n")

    data = pyobjcryst.powderpattern.PowderPattern()
    data.SetWavelength(wavelength)
    data.ImportPowderPattern2ThetaObs(pattern,0) # skip no lines
    ttheta = data.GetPowderPatternX()
    iobs = data.GetPowderPatternObs()

    # Create diffraction data, set crystal and input parameters for diffraction data
    diffData = data.AddPowderPatternDiffraction(crystal)
    diffData.SetReflectionProfilePar(pyobjcryst.powderpattern.ReflectionProfileType.PROFILE_PSEUDO_VOIGT,peakwidth * peakwidth)
    data.Prepare()
    icalc = data.GetPowderPatternCalc()

    # Check whether the largest peak of omitted 2theta range is at least as large as the largest peak of the observed 2theta range divided by PEAK_FACTOR
    warning=False
    full_ttheta,full_t1 = None,None
    if check_peaks and pattern is not None:
        PEAK_FACTOR = 3 # This is arbitrary,
        step = np.round((ttheta[1]-ttheta[0])/deg,3)
        full_ttheta = np.arange(0, max(ttheta/deg)+step, step)*deg
        data.SetPowderPatternX(full_ttheta)
        data.Prepare() # prepare for calcs, generate hkl space
        full_t1 = data.GetPowderPatternCalc()
        if max(full_t1[full_ttheta<min(ttheta)]) > max(icalc)/PEAK_FACTOR:
            warnings.warn("There is a significant diffraction peak in the omitted 2theta range! Statistical comparison can give biased results.", UserWarning)
            warning=True

    return ttheta,iobs,icalc,warning,full_ttheta,full_t1


def calculate_patterns(crystals, wavelength, peakwidth, pattern, check_peaks):
    # Basic argument check
    if (wavelength < 1.e-6):
        print("Warning: Wavelength is too small. A crash may be iminent.\n")

    data = pyobjcryst.powderpattern.PowderPattern()
    data.SetWavelength(wavelength)
    data.ImportPowderPattern2ThetaObs(pattern,0) # skip no lines
    ttheta = data.GetPowderPatternX()
    p = data.GetPowderPatternObs()

    patterns = np.zeros((len(crystals),len(ttheta)))
    for n,crystal in enumerate(crystals):
        data = pyobjcryst.powderpattern.PowderPattern()
        data.SetWavelength(wavelength)
        data.ImportPowderPattern2ThetaObs(pattern,0) # skip no lines
        ttheta = data.GetPowderPatternX()

        # Create diffraction data, set crystal and input parameters for diffraction data
        diffData = data.AddPowderPatternDiffraction(crystal)
        diffData.SetReflectionProfilePar(pyobjcryst.powderpattern.ReflectionProfileType.PROFILE_PSEUDO_VOIGT,peakwidth * peakwidth)
        data.Prepare()
        patterns[n] = data.GetPowderPatternCalc()

    return ttheta,p,patterns

def compare(pattern1, pattern2, plot, scale=True, warning=None, full_data=None):
    data = pyobjcryst.powderpattern.PowderPattern()
    ind = 0

    with open(pattern1,'r') as f:
        for line in f.readlines():
            ind += line.startswith('#')

    data.ImportPowderPattern2ThetaObs(pattern1,ind) # skip no lines
    ttheta1 = data.GetPowderPatternX()
    p1 = data.GetPowderPatternObs()

    ind = 0
    with open(pattern2,'r') as f:
        for line in f.readlines():
            ind += line.startswith('#')

    data.ImportPowderPattern2ThetaObs(pattern2,ind) # skip no lines
    ttheta2 = data.GetPowderPatternX()
    p2 = data.GetPowderPatternObs()

    # Calculate scale factor (p2 is reference)
    if scale:
        scalefactor = p2.max()/p1.max()
        scalefactor = FitScaleFactorForRw(p1,p2,scalefactor)
        p1 = p1 *scalefactor
        if full_data is not None:
            full_data[:,1] *= scalefactor

    # Check if xranges are equal
    try:
        assert (ttheta1==ttheta2).all() # this should always be True if the average pattern was derived using patterns outputted by this script
    except AssertionError:
        raise NotImplementedError('Functionality for non equal xranges has not yet been implemented!')
        # TODO: make interpolation function that can account for differing xranges in imported patterns

    # Compare data
    statistical_comparison(ttheta1,p1,p2,warning=warning)

    # Plot data
    if full_data is not None:
        full_data[:,0] *= deg

    if plot:
        plot_data(ttheta1,p1,p2,label1='pattern1',label2='pattern2', delta=False, warning=warning,full_data=full_data)


def background(pattern, peakwidth, bkg_points, locs, bkg_range, uniform, plot):
    assert bkg_points > 1
    data = pyobjcryst.powderpattern.PowderPattern()

    # Add the observed data
    data.ImportPowderPattern2ThetaObs(pattern,0) # skip no lines
    ttheta = data.GetPowderPatternX()

    # Background
    if locs is not None:
        bx = np.array(locs)*deg
    else:
        if bkg_range is not None:
            assert len(bkg_range)==2
            bx=np.linspace(bkg_range[0]*deg,bkg_range[1]*deg,bkg_points)
        else:
            bx=np.linspace(ttheta.min(),ttheta.max(),bkg_points)
        if not uniform:
            # adapt bx to minima of exp pattern in each neighbourhood (optional)
            ttheta = data.GetPowderPatternX()
            exp = data.GetPowderPatternObs()
            idx = [np.argmin(np.abs(ttheta-bxi)) for bxi in bx]

            step = (idx[1] - idx[0])//4
            for n in range(len(bx)):
                mn = -step if idx[n]>step else 0
                mx = step if n<(len(idx)-1) else 0
                bx[n] = ttheta[idx[n]+ mn + np.argmin(exp[idx[n]+mn:idx[n]+mx])]

    by=np.zeros(bx.shape)
    b=data.AddPowderPatternBackground()
    b.SetInterpPoints(bx,by)
    b.Print()
    b.UnFixAllPar()
    b.OptimizeBayesianBackground()

    no_bg = data.GetPowderPatternObs()-data.GetPowderPatternCalc()
    no_bg -= np.min(no_bg)

    if plot:
        plot_data(ttheta,data.GetPowderPatternCalc(),data.GetPowderPatternObs(),label1='background',label2='observed',vlines=bx,delta=True)

    with open('exp_no_bg.tsv','w') as f:
        for i in range(len(ttheta)):
            f.write("{:4.3f}\t{:10.8f}\n".format(ttheta[i]/deg,no_bg[i]))


def calculate(filename, crystal, wavelength, peakwidth, numpoints, max2theta, obspattern, plot, check_peaks, detail):
    data = pyobjcryst.powderpattern.PowderPattern()

    # Basic argument check
    if (wavelength < 1.e-6):
        print("Warning: Wavelength is too small. A crash may be iminent.\n")

    # Create diffraction data, set crystal and input parameters for diffraction data
    diffData = data.AddPowderPatternDiffraction(crystal)
    diffData.SetReflectionProfilePar(pyobjcryst.powderpattern.ReflectionProfileType.PROFILE_PSEUDO_VOIGT,peakwidth * peakwidth)

    # Set input parameters for pxrd
    data.SetWavelength(wavelength)

    # Add the observed data
    if obspattern is None:
        data.SetPowderPatternX(np.linspace(0, max2theta, numpoints))
        data.SetPowderPatternObs(np.ones(numpoints)) # Use fake unit intensity observed pattern since no observed pattern
    else:
        data.ImportPowderPattern2ThetaObs(obspattern,0) # skip no lines

    data.Prepare() # prepare for calcs, generate hkl space

    # Calculate initial scale factor
    if obspattern is not None:
        t1 = data.GetPowderPatternCalc()
        t2 = data.GetPowderPatternObs()
        scalefactor = t2.max()/t1.max()
        scalefactor = FitScaleFactorForRw(t1,t2,scalefactor)
    else:
        # don't scale data is there is no obspattern
        scalefactor = 1

    # Output data
    output_data(data,filename, detail)

    ttheta = data.GetPowderPatternX()
    iobs = data.GetPowderPatternObs()
    icalc = data.GetPowderPatternCalc()*scalefactor

    if obspattern is  None:
        # Plot data
        if plot:
            plot_data(ttheta,icalc,None)
    else:
        # Check whether the largest peak of omitted 2theta range is at least as large as the largest peak of the observed 2theta range divided by PEAK_FACTOR
        warning=False
        full_ttheta,full_t1 = None,None
        if check_peaks:
            PEAK_FACTOR = 3 # This is arbitrary,
            step = np.round((ttheta[1]-ttheta[0])/deg,3)
            full_ttheta = np.arange(0, max(ttheta/deg)+step, step)*deg
            data.SetPowderPatternX(full_ttheta)
            data.Prepare() # prepare for calcs, generate hkl space
            full_t1 = data.GetPowderPatternCalc()*scalefactor
            if max(full_t1[full_ttheta<min(ttheta)]) > max(icalc)/PEAK_FACTOR:
                warnings.warn("There is a significant diffraction peak in the omitted 2theta range! Statistical comparison can give biased results.", UserWarning)
                warning=True

        # Compare data - only sensical if there is a reference pattern
        statistical_comparison(ttheta,icalc,iobs,warning=warning)

        # Plot data
        if plot:
            plot_data(ttheta,icalc,iobs,warning=warning,full_data=np.array([full_ttheta,full_t1]).T)



def output_data(data, filename, detail):
    #Export data - calculated reflections
    output_name = filename.split('.')[0]
    calc = data.GetPowderPatternComponent(0)

    stdout_fd = sys.stdout.fileno()
    with open(output_name + '_fhkl.dat', 'w') as f, stdout_redirected(f):
        calc.PrintFhklCalc()

    if detail:
        with open(output_name + '_fhkl_detail.dat', 'w') as f, stdout_redirected(f):
            calc.PrintFhklCalcDetail()

    # Export data - 2theta space
    ttheta = data.GetPowderPatternX()
    iobs = data.GetPowderPatternObs()
    icalc = calc.GetPowderPatternCalc()*data.GetScaleFactor(0)

    with open(output_name + '.dat','w') as f:
        f.write('# 2theta \t ICalc \n')
        for i in range(len(ttheta)):
            f.write("{:4.3f}\t{:10.8f}\n".format(ttheta[i]/deg,icalc[i]))

    # Calculate - q space
    wavelength = calc.GetWavelength()
    qspace = 4.*np.pi/wavelength*np.sin(ttheta/2.)
    with open(output_name + '_q.dat','w') as f:
        f.write('# Q \t ICalc \n')
        for i in range(len(ttheta)):
            f.write("{:6.5f}\t{:10.8f}\n".format(qspace[i],icalc[i]))


def sp_output(ttheta,icalc,idx,warning,fttheta,ficalc):
    data = np.array([ttheta/deg, icalc]).T
    fdata = np.array([fttheta/deg, ficalc]).T

    if warning:
        data = np.vstack((data,[np.nan,np.nan]))

    np.save('frames/norm_{}.npy'.format(idx), data)
    np.save('frames/full_{}.npy'.format(idx), fdata)


##########################################

if __name__ == "__main__":
    print("This is a utils script, please run gpxrd.py!")
    sys.exit(1)
