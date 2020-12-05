#!/usr/bin/env python3.8
"""PSF convolution

  usage:
    convolv.py [-h|--help] -p psf1.fits -c param.json -o output.fits

  options:
   --help           show this help message and exit
   -p psf.fits
   -c param.json
   -o output.fits
"""
from docopt import docopt             # command line interface
import math
import numpy as np
import astropy.io.fits as fits
import json
import os
from scipy import ndimage 

#  command line interface
if __name__ == '__main__':
  args = docopt(__doc__)
  
  # convolution parameter
  with open(args['-c']) as f:
    p = json.load(f)
  
  ACEstd = int(p['ACEstd']['val']) # ACE [mas]
  
  # PSF
  hdul = fits.open(args['-p'])
  hdr = hdul[0].header
  psf = hdul[0].data
  
  M       = hdr['M']       # Number of FFT cells per wavelength in um
  
  cellscale = 1e-3 / M * 180*3600*1000 / math.pi # mas unit
  sigma = ACEstd / cellscale
  
  hdr['ACEstd'] = ACEstd
  hdr['ACEpix'] = sigma
  
  cpsf =  ndimage.gaussian_filter(psf,sigma=sigma)
  nhdu = fits.PrimaryHDU(cpsf)
  nhdu.header = hdr
  nhdulist = fits.HDUList([nhdu])
  nhdulist.writeto(args['-o'],overwrite=True)

