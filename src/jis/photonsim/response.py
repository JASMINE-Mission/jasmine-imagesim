import numpy as np
import sys
import math
from scipy import interpolate
from scipy import integrate
from scipy.constants import c, h
from jis.photonsim.fluxdensity import flux_rJHKs_byTeff, absmags


# wavelength definition of J and H bands (um).
WL_J = 1.250
WL_H = 1.644

def fluxdensity_spline(t_eff):
    """
    This function returns spline function for absolute magnitude
    instead of zero-mag flux of Vega-type star.

    Args: t_eff: effective temperature of star.

    Returns: spline function of flux density   
    """

    wav, flux   = flux_rJHKs_byTeff(t_eff)
    nphoton     = flux * wav * 1.0e-6 /h/c
    logL        = np.log10(wav)
    logN        = np.log10(nphoton)
    Npspline    = interpolate.interp1d(logL,logN,kind="cubic")

    return Npspline

def calc_response(control_params, telescope, detector):
    """
    This function calculates the electron rate (e-/s/m^2) detected by SJ
    based on the optics efficiency (EPdefined) and the quantum efficiency (QEdet).
    The target object is assumed to have SED determined by effective temperature t_eff.
    Hw-mag is calculated from the apparent J- and H-band magnitudes 
    using the same equation in the telescope_baseline as follows,
      Hw - H = 0.9*(J - H) - 0.06*(J - H),
    where J and H are the J- and H-band magnitudes in Vega system.

    Args:
        control_params: control parameters
        telescope: telescope object
        detector: detector object

    Returns:
        Tr  (float)  : Total electron rate (e-/s/m^2).
        WL  (ndarray): Wavelength grid of Npr (um; 0.1-um grid).
        Npr (ndarray): Electron flux (e-/s/m^2/um).

    Example:
       import numpy as np
       from jis.photonsim.response import calc_response

       JH  = 0.3
       alp = 0.75
       WL  = np.array([1.4, 1.5, 1.6])
       EP  = np.array([0.9, 0.8, 0.7])
       QEdet = np.array([0.5, 0.4, 0.3])

       el_rate, wavelength, el_flux = 


    """
    #Rv = control_params.Rv
    JH = control_params.JH
    #alp = control_params.alpha
    t_eff   = 9500
    #t_eff   = 3500
    #t_eff   = 6000
    WLdefined = telescope.opt_efficiency.wavelength
    EPdefined = telescope.opt_efficiency.efficiency
    WLshort = np.min(telescope.opt_efficiency.wavelength)
    WLlong  = np.max(telescope.opt_efficiency.wavelength)
    WLdet = detector.qe.wl
    QEdet = detector.qe.val

    _, J_abs, H_abs = absmags(t_eff)
    coeff1  = 0.9
    coeff2  = -0.06

    Ah  = Awl_n20(JH, J_abs, H_abs, WL_H)
    m_scale  = H_abs + (coeff1*JH + coeff2*JH**2 + Ah)

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
    EPinter = interpolate.interp1d(WLdefined, EPdefined, kind='linear',\
              bounds_error=False, fill_value=0.)
    QEinter = interpolate.interp1d(WLdet, QEdet, kind='linear',\
              bounds_error=False, fill_value=0.)
    EP = EPinter(WL)
    QE = QEinter(WL)

    # absolute-mag photon flux (ph/s/m^2/um).
    Np=Nphotons(WL, t_eff)

    # Photon count per wavelength in Hw-band.
    # The value is adjust by the distance-factor 
    # to make consistecy among the distance, extinction, and apparent Hw-mag.
    Npr = np.empty(len(WL))
    for i in range(len(WL)):
        Npr[i] = EP[i]*QE[i]*Np[i]*math.pow(10.0,-Awl_n20(JH,J_abs,H_abs,WL[i])/2.5)\
             * math.pow(10,m_scale/2.5)

    # Total photon (electron) rate (e-/s/m^2) from an object
    # with the same photon flux at the Hw band as a Zero-mag object.
    Tr = integrate.simps(Npr,WL)

    return Tr, WL, Npr 



def Awl_n20(JH, J_abs, H_abs, wl):
    '''
        Extinction A(wl) based on Nogueras-Lara et al.(2020;2020A&A...641A.141N).
        Yano-san derived an Awl relation based on NL20 as follows,
            Awl/A_K = 5.2106 * lambda^(-2.112).
        Also, JH - (M_J - M_H)  = E(J-H) = A_J - A_H
        A_K  = E(J-H)/(A_J/A_K - A_H/A_K)

        Awl = E(J-H)*wl^{-2.112}/(WL_J^-2.112 + WL_H^-2.2112)

    '''
    EJH = JH - (J_abs - H_abs) # E(J-H)
    c1  = 5.2106
    c2  = -2.112
    Ah_Ak   = c1*(WL_H**c2) # A_H / A_K
    Aj_Ak   = c1*(WL_J**c2) # A_J / A_K
    Ak  = EJH/(Aj_Ak - Ah_Ak)

    Awl = c1*(wl**c2)*Ak
    return Awl


def Nphotons(WL, t_eff):
    """
    This function returns logarithmically
    interpolated zero-mag photon flux (ph/s/m^2/um)
    in the MKO-NIR system at the wavelengths WL.

    Args:
        WL  (ndarray): Wavelength data (um).

    Returns:
        npw (ndarray): absolute-mag photon flux (ph/s/m^2/um).

    Example:
        import numpy as np
        from jis.photonsim.response import Nphotons

        wavelength = np.array([1.4, 1.5, 1.6])
        zm_ph_flux = Nphotons(wavelength, t_eff)

    """
    logwl = np.log10(WL)
    Npspline  = fluxdensity_spline(t_eff)
    lognp = Npspline(logwl)    # 
    npw   = np.power(10,lognp)

    return npw
