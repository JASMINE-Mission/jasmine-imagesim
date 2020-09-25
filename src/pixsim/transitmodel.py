import numpy as np
import argparse
import batman
import astropy.units as u
from astropy.constants import G, R_sun, M_sun, R_earth, M_earth
import numpy as np
import matplotlib.pyplot as plt

def inject_transit(p,t0,depth, ntime):
    params = batman.TransitParams()
    params.t0 = t0                      #time of inferior conjunction
    params.per = p                      #orbital period
    params.rp = np.sqrt(depth)                      #planet radius (in units of stellar radii)
    params.a = 15.                       #semi-major axis (in units of stellar radii)
    params.inc = 87.                     #orbital inclination (in degrees)
    params.ecc = 0.                      #eccentricity
    params.w = 90.                       #longitude of periastron (in degrees)
    params.u = [0.1, 0.3]                #limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"       #limb darkening model
    
    t = np.linspace(0, ntime, ntime)
    m = batman.TransitModel(params, t)    #initializes model
    flux = m.light_curve(params)          #calculates light curve
    return flux

def gentransit(t,t0=0.0,Porb=14.0,Rp=1.0,Mp=1.0,Rs=0.2,Ms=0.15,ideg=90.0,w=90.0,e=0.0,u1=0.1,u2=0.3):
    #Rp Earth radius
    params = batman.TransitParams()
    params.t0 = t0 # time of inferior conjunction
    params.rp = Rp*R_earth/(Rs*R_sun) # planet radius (in units of stellar radii)
    
    # calculate semi-major axis from orbital period value                        
    params.inc = ideg  # orbital inclination (in degrees)                        
    params.ecc = e # eccentricity                                                
    params.w = w # longitude of periastron (in degrees), 90 for circular          
    params.u = [u1,u2] # limb darkening coefficients                              
    params.limb_dark = "quadratic" # limb darkening model                        
    
    #period update                                                                
    params.per = Porb # orbital period (days)                                    
    a = (((params.per*u.day)**2 * G * (Ms*M_sun + Mp*M_earth) / (4*np.pi**2))**(1./3)).to(R_sun).value/Rs
    params.a = a # semi-major axis (in units of stellar radii)                    
    b=a*np.cos(ideg/180.0*np.pi)
    m = batman.TransitModel(params, t) # initializes the model                    
    injlc = np.array(m.light_curve(params))

    return injlc, b

