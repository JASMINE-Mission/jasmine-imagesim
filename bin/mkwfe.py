#!/usr/bin/env python3.8

""" Make a wave front error map.

.. code-block:: bash

  usage:
    mkwfe.py [-h|--help] -t tel.json -e wfe.json -m wfe.fits

  options:
    -h --help    show this help message and exit
    -t tel.json  telescope parameter file
    -e wfe.json  wave front error parameter file
    -m wfe.fits  output wave front error map
"""
from docopt import docopt             # command line interface
import astropy.io.fits as fits
import json
import sys
from jis.photonsim import extract_json
from jis.photonsim import wfe
from jis.photonsim import readfits

#  Command line interface
if __name__ == '__main__':
  args = docopt(__doc__)
  
  # Get telescope parameter
  with open(args['-t']) as f:
    q = json.load(f)
    EPD = float(q['EPD']['val'])  
  
  data=wfe.calc_wfe(EPD,args['-e'])
  
  # Save WFE map
  hdu = fits.PrimaryHDU(data)
  hdu.header["WFE-FILE"] = args['-e']
  hdu.header["WFE-EPD"] = EPD
  hdulist = fits.HDUList([hdu])
  hdulist.writeto(args['-m'],overwrite=True)
