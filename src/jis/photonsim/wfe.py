import math
import numpy as np
from jis import photonsim
from jis.photonsim import zernike                       # Zernike polynomials functions
import json

def calc_wfe(EPD,efile):
    """
    Summary:
        This function calculates wavefront error pattern 
        by calculating summation of zernike polynomials.
        Coefficients and dimensions of the zernike terms are 
        given by the input json file, efile.

    Args:
        EPD   (float): Entrance pupil diameter (pix=mm).
        efile (str)  : Filename of the json file having wavefront error parameters.
        *** Pixel scale is assumed to be 1 mm/pix!!! ***

    Returns:
        data (ndarray): EPD+4 x EPD+4 data array of the calculated wavefront error pattern.

    Example:
        from jis.photonsim.wfe import calc_wfe

        wfe_pattern = calc_wfe(500, "wfe.json").

    """

    # make a little bit larger data
    N    = int(EPD + 4) 
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
