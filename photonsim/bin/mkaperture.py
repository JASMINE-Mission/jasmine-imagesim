#!/usr/bin/env python3.8
"""
  Make an aperture image

  usage:
    mkaperture.py [-h|--help] -t tel.json -a aperture.fits

 options:
   -h --help        show this help message and exit
   -t tel.json      telescope parameter file
   -a aperture.fits output aperture image
"""
from docopt import docopt             # command line interface
import math
import numpy as np
import astropy.io.fits as fits
import json
import sys

#  command line interface
if __name__ == '__main__':
  args = docopt(__doc__)

# get telescope parameter
with open(args['-t']) as f:
  p = json.load(f)

EPD = float(p['EPD']['val'])  # D must be given in mm unit. A cell size is 1 mm.
Cobs = float(p['Cobs']['val'])
Robs = EPD/2*Cobs
Stype = p['Stype']['val']
if Stype  == 'tripod':
  Tsp =  float(p['Stype']['thick'])
else:
  print('Stype error')
  sys.exit()

N = int(EPD+4)   # 2mm larger than D
if N%2 == 1:   # should be even
  N = N+1 

# initialize data
data=np.zeros((N,N),dtype=np.float64) # 

# Set aperture
a = math.sqrt(3)/2
S=0
for i in range(N):
  y = i - N/2 # y axis,  origin at N/2
  for j in range(N):
    x = j - N/2   # x axis
    r = x*x + y*y # square of distance from the center 
    x1 = -x/2 + y*a  # x1,y1 is rotated -120 degree 
    y1 = -x*a - y/2
    x2 = -x/2 - y*a  # x2,y2 is rotated 120 degree 
    y2 =  x*a - y/2
    apt = 0.0 
    if r <= (EPD/2.)**2  and r > Robs**2 : 
      apt = 1.0 
      if x  >0 and y  > -Tsp/2 and y  < Tsp/2: # in a spider
        apt = 0.0
      if x1 >0 and y1 > -Tsp/2 and y1 < Tsp/2: # in a spider
        apt = 0.0
      if x2 >0 and y2 > -Tsp/2 and y2 < Tsp/2: # in a spider
        apt = 0.0
    #NOTE: if Rob==0 and Tsp==0 then apt=1 for r==0 and y==0 
    data[i,j] = apt 
    S = S + apt

# save the data
hdu = fits.PrimaryHDU(data)
hdu.header['APTFILE'] = args['-t']  # telescope parameter file
hdu.header['EPD']      = EPD         # pupil diameter (mm)
hdu.header['COBS']     = Cobs        # central obscuration ratio 
hdu.header['STYPE']    = Stype       # Spider type
hdu.header['TSP']      = Tsp         # spider thickness (mm)
hdu.header['STEL']     = S*1e-6      # Area of telescope aperture (m2)
hdulist = fits.HDUList([hdu])
hdulist.writeto(args['-a'],overwrite=True)
