#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""Make a PSF

.. code-block:: bash

  usage:
    mkpsf.py [-h|--help] -a apt -w wfe -s spec -c cntl -n fn -p psf.fits

  options:
    -h --help    show this help message and exit
    -a apt       aperture mask
    -w wfe       wfe map or null
    -s spec      electrons/sec/um at lambada
    -c cntl      control parameter
    -n fn        number of fp-cells of the output psf image.
    -p psf.fits
"""
from docopt import docopt             # command line interface
import json
import sys
from jis import photonsim
from jis.photonsim import extract_json
from jis.photonsim import psf
from jis.photonsim import readfits
import astropy.io.fits as fits

#  command line interface
if __name__ == '__main__':
  args = docopt(__doc__)
  
  # spectral response
  with open(args['-s']) as f:
    sp = json.load(f)
  
  k,WL,NP,Ntot=extract_json.extsp(sp)
  Stel,adata,aN,ahdr=readfits.read_aperture_mask(args['-a'])
  
  # get the WFE map
  if args['-w'] == 'null':
    wN = aN
    wfe = np.zeros((wN,wN))
  else:
    wfe,wN,whdr=readfits.read_wfe_map(args['-w'])
  
  # control parameter
  with open(args['-c']) as f:
    cp = json.load(f)
  M = cp['M']['val']
 
  # get the number of fp-cells
  if args['-f']:
    fN = int(args['-f'])
  else:
    fN = None

  ### psf
  if fN is None:
    image = psf.calc_psf(wfe, wN, k, WL, NP, Ntot, Stel, adata, M, aN)
  else:
    image = psf.calc_psf(wfe, wN, k, WL, NP, Ntot, Stel, adata, M, aN, fN)
  
  # Save 
  hdu = fits.PrimaryHDU(image)
  hdu.header['NTOT'] = Ntot  # total flux/m2/sec
  
  hdu.header['APTFILE'] = ahdr['APTFILE']     # telescope.json
  hdu.header['EPD']     = ahdr['EPD']         # pupil diameter (mm)
  hdu.header['COBS']    = ahdr['COBS']        # central obscuration ratio 
  hdu.header['STYPE']   = ahdr['STYPE']       # Spider type
  hdu.header['TSP']     = ahdr['TSP']         # spider thickness (mm)
  hdu.header['STEL']    = ahdr['STEL']        # Area of telescope aperture (m2)
  hdu.header['SPRES']   = args['-s']          # spectral responce 
  hdu.header['M'] = M   # Number of FFT cells per wavelength in um
  
  hdulist = fits.HDUList([hdu])
  hdulist.writeto(args['-p'],overwrite=True)
  
