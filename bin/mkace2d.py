#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""Make an atitude control error map

  usage:
    mkace2d.py [-h|--help] -x aX.fits -y aY.fits -v xs -w ys -p ps -n N -o out.fits

 options:
   --help       show this help message and exit
   -x aX.fits   X-axis simulated ACE file
   -v xs        aX.fits is scaled by xs
   -y aY.fits   Y-axis simulated ACE file
   -w ys        aY.fits is scaled by ys
   -n N         output image is N x N 
   -p ps        pixel scale of Output image
   -o out.fits  output image file

""" 
from docopt import docopt             # command line interface
import math
import numpy as np
import astropy.io.fits as fits
import sys
from jis.photonsim import ace

# Command line interface
if __name__ == '__main__':
  args = docopt(__doc__)
  
  # Get parameters from command line
  x_scale = float(args['-v'])
  y_scale = float(args['-w'])
  pix_scale = float(args['-p'])
  N = int(args['-n'])
  
  xhdul = fits.open(args['-x'])
  xdata = xhdul[0].data
  xhead = xhdul[0].header
  xN = xhead['NAXIS1']
  yhdul = fits.open(args['-y'])
  ydata = yhdul[0].data
  yhead = yhdul[0].header
  yN = xhead['NAXIS1']
  
  if xN != yN:
    print("NAXIS1 is not the same")
    sys.exit(-1)
  
  data = ace.ace2d(x_scale,y_scale,pix_scale,N,xdata,ydata,yN)
  
  hdu = fits.PrimaryHDU(data)
  hdulist = fits.HDUList([hdu])
  hdulist.writeto(args['-o'],overwrite=True)
