import astropy.io.fits as fits
import sys
import math

def read_wfe_map(fitsfile):
  whdul = fits.open(fitsfile)
  whdr = whdul[0].header
  wN = whdr['NAXIS1']
  if whdr['NAXIS'] != 2:
    #    print("NAXIS = {0} is bad in {1}".format(whdr['NAXIS'],args['-e']))
    print("NAXIS = {0} is bad in {1}".format(whdr['NAXIS']))
    sys.exit(-1)
  if whdr['NAXIS1'] != whdr['NAXIS2']:
    print("NAXIS1 {0} != NAXIS2 {1}".format(wN,whdr['NAXIS2']))
    sys.exit(-1)
  wfe = whdul[0].data *2*math.pi
  return wfe,wN,whdr

def read_aperture_mask(fitsfile):
    # get the aperture mask
    ahdul = fits.open(fitsfile)
    ahdr = ahdul[0].header
    aN = ahdr['NAXIS1']
    if ahdr['NAXIS1'] != ahdr['NAXIS2']:
        print("NAXIS1 {0} != NAXIS2 {1}".format(aN,ahdr['NAXIS2']))
        sys.exit()
    Stel = ahdr['STEL']   # Area of collecting photons
    adata = ahdul[0].data
    return Stel,adata,aN,ahdr
