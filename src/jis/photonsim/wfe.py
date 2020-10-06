import math
import numpy as np
from jis import photonsim
from jis.photonsim import zernike                       # Zernike polynomials functions
import json

def calc_wfe(EPD,efile):
    # make a little bit larger data
    N = int(EPD + 4) 
    data = np.zeros((N,N),dtype=float)
    
    # Get wavefront error parameters from json file
    with open(efile) as f:
        p = json.load(f)
        NP = int(p['N-polys'])
        Zn = np.empty(NP)
        Zm = np.empty(NP)
        Za = np.empty(NP)
        for i in range(NP):
            Zn[i] = int(p['z{:03d}-n'.format(i+1)])
            Zm[i] = int(p['z{:03d}-m'.format(i+1)])
            Za[i] = float(p['z{:03d}-a'.format(i+1)])
            
            
    for i in range(NP):
        for iy in range(N):
            for ix in range(N):
                y = iy - N/2
                x = ix - N/2
                rho = math.sqrt( y*y+x*x ) / (EPD/2) 
                if rho <= 1 :
                    th = math.atan2( y, x)
                    data[iy,ix] = data[iy,ix] + Za[i]*zernike.Zernike(Zn[i],Zm[i],rho,th)
                else :
                    data[iy,ix] = math.nan

    return data
