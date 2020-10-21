import numpy as np
import argparse
import batman
import astropy.units as u
from astropy.constants import G, R_sun, M_sun, R_earth, M_earth
import numpy as np
import matplotlib.pyplot as plt

def inject_transit(p, t0, depth, ntime,\
                   a=15., inc=87., ecc=0., w=90.,\
                   ld_coeff=[0.1, 0.3], ld_model='quadratic'):
    """
    Summary:
        This function calculates a transit light curve.
        The transit model is calculated based on the input parameters.

    Args:
        p        (float): Orbital period.
        t0       (float): Time of inferior conjunction (same unit as p).
        depth    (float): Transit depth (area ratio of the planet and the central star).
        ntime    (int)  : Time of the end of the calculation (same unit as p).
        a        (float): Semi-major axis in Rstar.
                          (Default: 15.)
        inc      (float): Orbital inclination in deg.
                          (Default: 87.).
        ecc      (float): Eccentricity in deg.
                          (Default: 0.).
        w        (float): Longitude of periastron in deg.
                          (Default: 90.)
        ld_coeff (array): Limb-darkening coefficients [u1, u2]
                          (Default: [0.1, 0.3])
        ld_model (str)  : Limb-darkening model
                          (Default: quadratic).

    Returns:
        flux (ndarray): Transit light curve data.

    """

    params = batman.TransitParams()
    params.t0  = t0                       # time of inferior conjunction
    params.per = p                        # orbital period
    params.rp  = np.sqrt(depth)           # planet radius (in units of stellar radii)
    params.a   = a                        # semi-major axis (in units of stellar radii)
    params.inc = inc                      # orbital inclination (in degrees)
    params.ecc = ecc                      # eccentricity
    params.w   = w                        # longitude of periastron (in degrees)
    params.u   = ld_coeff                 # limb darkening coefficients [u1, u2]
    params.limb_dark = ld_model           # limb darkening model
    
    t = np.linspace(0, ntime, ntime)
    m = batman.TransitModel(params, t)    # initializes model
    flux = m.light_curve(params)          # calculates light curve

    return flux


def gentransit(t, t0=0.0, Porb=14.0, Rp=1.0, Mp=1.0, Rs=0.2, Ms=0.15,\
               ideg=90.0,w=90.0,e=0.0,u1=0.1,u2=0.3):
    """
    Summary:
        This function calculates a transit light curve
        based on the input parameters.

    Args:
        t      (ndarray): Array of times at which 
                          to calculate the transit model in days.
        t0     (float)  : Time of inverior conjunction in days.
                          (Default: 0.0)
        Porb   (float)  : Orbital period in days.
                          (Default: 14.0)
        Rp     (float)  : Planet radius in Earth radius.
                          (Default: 1.0)
        Mp     (float)  : Planet mass in Earth mass.
                          (Default: 1.0)
        Rs     (float)  : Stellar radius in Solar radius.
                          (Default: 0.2)
        Ms     (float)  : Stellar mass in Solar mass.
                          (Default: 0.15)
        ideg   (float)  : Orbital inclination in deg.
                          (Default: 90.0)
        w      (float)  : Longitude of periastron in deg.
                          (Default: 90.0)
        e      (float)  : Eccentricity in deg.
                          (Default: 0.0)
        u1, u2 (float)  : Limb-darkening parameters.
                          (Default: u1=0.1; u2=0.3)

    Returns:
        injlc (ndarray): Transit light curve data.
        b     (float)  : Apparent size of the orbital semi-major axis in Rstar? 
          
    """

    params = batman.TransitParams()

    params.t0 = t0 # time of inferior conjunction

    # Rp in Earth radius.
    # Convert the unit and set it to rp in params.
    params.rp = Rp*R_earth/(Rs*R_sun) # planet radius (in units of stellar radii)
   
    # Set orbital parameters and limb-darkening parameters. 
    params.inc = ideg    # orbital inclination (in degrees)                        
    params.ecc = e       # eccentricity                                                
    params.w   = w       # longitude of periastron (in degrees), 90 for circular          
    params.u   = [u1,u2] # limb darkening coefficients                              
    params.limb_dark = "quadratic" # limb darkening model                        
    
    # period update                                                                
    params.per = Porb # orbital period (days)                                    

    # calculate semi-major axis from the orbital period parameters.
    a = (((params.per*u.day)**2 * G * (Ms*M_sun + Mp*M_earth) / (4*np.pi**2))**(1./3)).to(R_sun).value/Rs
    params.a = a # semi-major axis (in units of stellar radii)                    

    b=a*np.cos(ideg/180.0*np.pi)

    m     = batman.TransitModel(params, t)  # initializes the model                    
    injlc = np.array(m.light_curve(params)) # get the light curve data

    return injlc, b

