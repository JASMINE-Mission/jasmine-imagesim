#!/usr/bin/env python3.8

"""
  Make a wave fromt error map

  usage:
    mkwfe.py [-h|--help] -t tel.json -e wfe.json -m wfe.fits

 options:
   --help       show this help message and exit
   -t tel.json  telescope parameter file
   -e wfe.json  wave front error parameter file
   -m wfe.fits  output wave front error map
"""
from docopt import docopt             # command line interface
import math
import numpy as np
import astropy.io.fits as fits
import json
import sys
import zernike                       # Zernike polynomials functions

#  Command line interface
if __name__ == '__main__':
  args = docopt(__doc__)

# Get telescope parameter
with open(args['-t']) as f:
  q = json.load(f)
  EPD = float(q['EPD']['val'])  

# make a little bit larger data
N = int(EPD + 4) 
data = np.zeros((N,N),dtype=float)

# Get wavefront error parameters from json file
with open(args['-e']) as f:
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

# Save WFE map
hdu = fits.PrimaryHDU(data)
hdu.header["WFE-FILE"] = args['-e']
hdu.header["WFE-EPD"] = EPD
hdulist = fits.HDUList([hdu])
hdulist.writeto(args['-m'],overwrite=True)
