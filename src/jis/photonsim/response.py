import numpy as np
import sys
import math
from scipy import interpolate
from scipy import integrate   

# physical constants
c = 299792458.     ; # [m/s] 
h = 6.62607015e-34 ; # [J/s] 

# MKO NIR system from 2005PASP..117..421T, 2005PASP..117.1459T
## wavelength [um]
WL_MKO   = np.array([ 0.5450 , 1.250  , 1.644  , 2.121  , 2.149  , 2.198  ,
                      3.754  , 4.702])
## zero-mag flux [W/m^2/um]
FLUX_MKO = np.array([3.68e-08,3.01e-09,1.18e-09,4.57e-10,4.35e-10,4.00e-10,
                     5.31e-11,2.22e-11])

# zero-mag photon flux [Photons/s/m^2/um]
NP_MKO = FLUX_MKO * WL_MKO * 1.0e-6 / h / c ;

# interpolate photon flux (photons/s/m^2/um)
logL = np.log10(WL_MKO)
logN = np.log10(NP_MKO)
Npspline = interpolate.interp1d(logL, logN, kind="cubic")

# A(lambda)/Av
AWL = lambda wl , R : EXT_a(wl) + EXT_b(wl)/R

# wavelength definition of J and H bands (um).
WL_J = 1.250
WL_H = 1.644

def calc_response(Rv, JH, alp, k, WLdefined, EPdefined, WLshort, WLlong, WLdet, QEdet):
    """
    This function calculates the electron rate (e-/s/m^2) detected by SJ
    based on the optics efficiency (EPdefined), the quantum efficiency (QEdet).
    The target object is assumed to have an interstellar extinction
    defined with Rv and JH (=E(J-H)) and show the same photon flux 
    at the Hw-band wavelength as a zero-mag object.
    The parameter 'alp' defines the weight to determine the Hw-band magnitude
    by interpolating the J- and H-band magnitudes.

    Args:
        Rv        (float)  : Extinction parasmeter Rv(=Av/E(B-V)).
        JH        (float)  : Color excess E(J-H)(=AJ-AH).
        alp       (float)  : Interpolation factor to define Hw-band mag.
                             Hw-band mag = alp * J-mag + (1-alp) H-mag
        k         (int?)   : Number of wavelength data points? (not used).
        WLdefined (ndarray): Wavelength data.
        EPdefined (ndarray): Optical efficiency data.
        WLshort   (float)  : Shortest wavelength (um).
        WLlong    (float)  : Longest wavelength (um). 
        WLdet     (ndarray): Wavelength grid for the detector QE data, QEdet.
        QEdet     (ndarray): Quantum efficiency data of the detector.

    Returns:
        Tr  (float)  : Total electron rate (e-/s/m^2).
        WL  (ndarray): Wavelength grid of Npr (um).
        Npr (ndarray): Electron flux (e-/s/m^2/um).

    Example:
       import numpy as np
       from jis.photonsim.response import calc_response

       Rv  = 3.1
       EJH = 0.3
       alp = 0.75
       WL  = np.array([1.4, 1.5, 1.6])
       EP  = np.array([0.9, 0.8, 0.7])
       QEdet = np.array([0.5, 0.4, 0.3])
       
       el_rate, wavelength, el_flux = \
           calc_response(Rv, EJH, alp, len(WL), WL, EP, np.min(WL), np.max(WL), WL, QEdet) 

    """
      
    # from Rv and J-H, calculate Av
    JHA = AWL(WL_J,Rv) - AWL(WL_H,Rv) # (AJ-AH)/Av
    Av  = JH/JHA
    
    # Array for wavelength
    # WL[i] should be in (WLshort,WLlong and 0.1 um step)
    Wlist=[]
    eps = 1e-8 
    for i in range(22):
        w = i/10+0.4  # from 0.4 um to 2.5 um
        if w >= WLshort-eps and w <= WLlong+eps :
            Wlist.append(w)
    WL = np.array(Wlist)

    # interpolate efficiency and qe
    EPinter = interpolate.interp1d(WLdefined, EPdefined, kind='linear')
    QEinter = interpolate.interp1d(WLdet, QEdet, kind='linear')
    EP = EPinter(WL)
    QE = QEinter(WL)

    # Zero-mag photon flux (ph/s/m^2/um)
    Np=Nphotons(WL)

    # band definition
    NpJ   = Nphotons(WL_J)      # Zero-mag photon flux in J  band (ph/s/m^2/um).
    NpH   = Nphotons(WL_H)      # Zero-mag photon flux in H  band (ph/s/m^2/um).
    NpHw  = NpJ*alp+NpH*(1-alp) # Zero-mag photon flux in Hw band (ph/s/m^2/um).
    ## Photon fluxes of reddened object which is intrinsically zero mag. 
    NprJ  = NpJ*math.pow(10.0,-AWL(WL_J,Rv)*Av/2.5)     # J-band reddened photon flux (ph/s/m^2/um).
    NprH  = NpH*math.pow(10.0,-AWL(WL_H,Rv)*Av/2.5)     # H-band reddened photon flux (ph/s/m^2/um).
    NprHw = NprJ*alp+NprH*(1-alp) # Hw-band reddened photon flux (ph/s/m^2/um).

    # reddenend photon (electron) flux (e-/s/m^2/um) from an object
    # with the same photon flux at the Hw band as a Zero-mag object.
    Npr = np.empty(len(WL))
    for i in range(len(WL)):
        Npr[i] = EP[i]*QE[i]*Np[i]*math.pow(10.0,-AWL(WL[i],Rv)*Av/2.5)*NpHw/NprHw

    # Total photon (electron) rate (e-/s/m^2) from an object
    # with the same photon flux at the Hw band as a Zero-mag object.
    Tr = integrate.simps(Npr,WL)

    return Tr, WL, Npr

def Nphotons(WL):
    """
    This function returns logarithmically 
    interpolated zero-mag photon flux (ph/s/m^2/um)
    in the MKO-NIR system at the wavelengths WL.

    Args:
        WL  (ndarray): Wavelength data (um).

    Returns:
        npw (ndarray): Zero-mag photon flux (ph/s/m^2/um).

    Example:
        import numpy as np
        from jis.photonsim.response import Nphotons

        wavelength = np.array([1.4, 1.5, 1.6])
        zm_ph_flux = Nphotons(wavelength) 

    """ 
    logwl = np.log10(WL)
    lognp = Npspline(logwl)    # globally defined in this file.
    npw   = np.power(10,lognp)

    return npw

def EXT_a(wl):
    """
    This function returns the a(x) coefficient defining the 
    interstellar extinction curve in 1989ApJ...345..245C.
    A(lmd)/Av = a(x) + b(x)/Rv; x = 1/lmd(um); 

    Args:
        wl (ndarray): Wavelength (um).

    Returns:
        ret (ndarray): Coefficients a(x).

    Example:
        import numpy as np
        from jis.photonsim.response import EXT_a

        wl = np.array([1.4, 1.5, 1.6])
        a_coeff = EXT_a(wl)

    """ 

    x = 1/wl

    if x >= 0.3 and x<=1.1 :
        ret = 0.574 * math.pow(x,1.61)
    elif x <= 8. :
        y = x - 1.82
        y1 = 1+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4
        y2 =   0.01979*y**5-0.77530*y**6+0.32999*y**7
        ret = y1 + y2 
    else :
        print('WL out of range error')
        sys.exit()

    return ret
    
def EXT_b(wl):
    """
    This function returns the b(x) coefficient defining the
    interstellar extinction curve in 1989ApJ...345..245C.
    A(lmd)/Av = a(x) + b(x)/Rv; x = 1/lmd(um);

    Args:
        wl (ndarray): Wavelength (um).

    Returns:
        ret (ndarray): Coefficients b(x).

    Example:
        import numpy as np
        from jis.photonsim.response import EXT_b

        wl = np.array([1.4, 1.5, 1.6])
        b_coeff = EXT_b(wl)

    """

    x = 1/wl

    if x >= 0.3 and x<=1.1:
        ret = -0.527* math.pow(x,1.61)
    elif x <= 8.:
        y  = x - 1.82
        y1 = 1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5
        y2 = 5.30260*y**6-2.09002*y**7
        ret = y1 + y2
    else :
        print('WL out of range error')
        sys.exit()

    return ret

