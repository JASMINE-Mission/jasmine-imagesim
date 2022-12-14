#!/usr/bin/env python3
import numpy as np
import sys
from matplotlib import pyplot as plt
import pkgutil
from io import BytesIO

from scipy import interpolate
from scipy.constants import c,h

# MKO NIR system from Tokunaga & Vacca 05 (2005PASP..117..421T, 2005PASP..117.1459T)
# wavelength [um] (isophotal wavelength for Vega model)
WL_MKO   = np.array([0.5450, 1.250, 1.644, 2.121, 2.149]) # V/J/H/K'/Ks
                      #, 2.198, 3.754, 4.702]) # K/L'/M'
# flux [W/m2/um] (mean flux/flux at the iso. wavelength of the Vega model)
FLUX_MKO = np.array([3.68e-08, 3.01e-09, 1.18e-09, 4.57e-10, 4.35e-10]) # V/J/H/K'/Ks
                      #, 4.00e-10, 5.31e-11, 2.22e-11]) # K/L'/M'

# Settings for the magnitudes on maglist.txt
# wavelength of the r, J, H, and Ks bands [um]
WL_rJHK  = np.array([0.6165, 1.250, 1.644, 2.149])

def ABmag_to_Flam(ABmag, l_p):
    """Transform AB magnitude to F_lamda

    Args:
        ABmag (float)   : AB magnitude
        l_p (float)     : pivot wavelngth (um) (see Tokunaga & Vacca, 2005)

    Returns:
        f_lam (float)   : Flux densities (W/m2/um)
    """

    f_nu    = 3720*10**(-0.4*ABmag) # Jy
    f_nu    = f_nu*1.0e-26          # W/m2/Hz
    f_lam   = f_nu/((l_p**2)/c*1e-6)

    return f_lam

def VEGAmag_to_ABmag(VEGAmag, Fo):
    """Transform Vega magnitude to AB magnitude

    Args:
        VEGAmag (float) : Vega magnitude
        Fo (float)      : photometric zero flux for Vega mag.

    Returns:
        ABmag (float)   : AB magnitude
    """
    F_nu    = Fo*10**(-VEGAmag/2.5) # Jy
    ABmag   = -2.5*np.log10(F_nu) + 8.926

    return ABmag

def F_rJHK():
    """Calculate F_lambda and ABmag in the r, J, H, and Ks bands
       for dwarf stars at 10 pc.

    Returns:
        f_lam (ndarray) : F_lambda (W/m2/um)
        t_eff (ndarray) : Effective temperature (K)
        ABmags (ndarray): AB magnitudes

    This caluculation is based on Kraus & Hillenbrand (2007)
    and Tokunaga & Bacca (2005).
    Note the SDSS photometry (here r-band) is intended to be
    on the AB system (Oke & Gunn 1983).
    lam_pvt_* is pivot wavelength descrived in Eq(A10) in
    Tokunaga & Vacca (2005).

    """
    lam_pvt_r   = 0.618

    #d_tab5  = np.loadtxt("photonsim/data/maglist.txt", comments='#', dtype='f8')
    maglist = pkgutil.get_data("jis", "photonsim/data/maglist.txt")
    d_tab5  = np.loadtxt(BytesIO(maglist), comments='#', dtype='f8')

    t_eff   = d_tab5[:,-1]
    mags_ar = d_tab5[:,:-1]

    r_pzero = 3631 #zero flux in r band [Jy]
    fo_ar   = (r_pzero, 1594, 1024, 666.7)
    ABmags  = VEGAmag_to_ABmag(mags_ar, fo_ar)

    # lambda_pivot for each bandpass (Table 2 in Tokunaga & Vacca 2005)
    lp_ar   = np.array((lam_pvt_r, 1.247, 1.628, 2.150), dtype='f8')
    f_lam   = ABmag_to_Flam(ABmags, lp_ar)

    return f_lam, t_eff, ABmags


def flux_rJHKs_byTeff(t_eff=5500):
    """Calculate F_lamda array in the r, J, H, and Ks bands
       for a dwarf star at 10 pc with a given effective temperature.

    Args:
        Teff (float)    : Effective temperature of target

    Returns:
        WL_rJHK         : Wavelengths of the r, J, H, and Ks bands
        FLUX_rJHK       : Flux densities in the r, J, H, and Ks bands (W/m2/um)

    WL_rJHK and FLUX_rJHK can be alternatives to the WL_MKO and FLUX_MKO
    which are currently used in response.py, respectively.
    These flux densities let the imageim test different spectral types
    for targets.
    """

    Flam, t_eff_ar, ABmags  = F_rJHK()
    Flam        = Flam.T
    ABmags      = ABmags.T

    f   = interpolate.interp1d(t_eff_ar, Flam, kind='cubic')
    FLUX_rJHK   = f(t_eff)

    return  WL_rJHK, FLUX_rJHK

def absmags(t_eff=5500):
    """
    derive absolute VEGA magnitude from T_eff.

    args:
        t_eff (float)   : effective temperature
    returns:
        hw_mag (float)  : absolute vega magnitude in Hw band
    """
    #d_tab5  = np.loadtxt("photonsim/data/maglist.txt", comments='#', dtype='f8')
    maglist = pkgutil.get_data("jis", "photonsim/data/maglist.txt")
    d_tab5  = np.loadtxt(BytesIO(maglist), comments='#', dtype='f8')
    teff_ar = d_tab5[:,-1]
    mags_ar = d_tab5[:,:-1].T

    f = interpolate.interp1d(teff_ar, mags_ar, kind='cubic')
    mag_ar  = f(t_eff)
    J_mag   = mag_ar[1]
    H_mag   = mag_ar[2]
    Hw_mag  = 0.9*J_mag + 0.1*H_mag - 0.06*(J_mag - H_mag)**2

    return Hw_mag, J_mag, H_mag


if __name__=='__main__':

    if len(sys.argv) <= 2:
        if len(sys.argv) == 2:
            t_eff   = float(sys.argv[1])
        else:
            t_eff   = 5500

        Hw_mag,_,_  = absmags(t_eff)
        print("absolute Hw mag:\t",Hw_mag)

        WL_rJHK, FL_rJHK     = flux_rJHKs_byTeff(t_eff)
        print(WL_rJHK, FL_rJHK)
        plt.plot(WL_rJHK, FL_rJHK)
        plt.savefig("fluxdensity_"+str(t_eff)+".png")
    else:
        print("args) [effective temperature]")
