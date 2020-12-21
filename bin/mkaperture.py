#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""Make an aperture image

.. code-block:: bash

  usage:
    mkaperture.py [-h|--help] -t tel.json -a aperture.fits

  options:
    -h --help         show this help message and exit
    -t tel.json       telescope parameter file
    -a aperture.fits  output aperture image
"""
from docopt import docopt             # command line interface
import math
import numpy as np
import astropy.io.fits as fits
import json
import sys
from jis import photonsim
from jis.photonsim import aperture

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
  
  data,S=aperture.calc_aperture(N,EPD,Robs,Tsp)
  
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
