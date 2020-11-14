#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""
  Make an image

  usage:
    mkimage.py [-h|--help] [--pd paramdir] --starplate star_plate.csv --det det.json --tel tel.json --ace ace.json --ctl ctl.json [--od outdir] [--overwrite] 

 options:
   --help                     show this help message and exit.
   --pd paramdir              name of the directory containing parameter files.
   --starplate star_plate.csv csv file containing star info (plate_id, star_id, xpix, ypix, l, b)
   --det det.json             json file containing detector related parameters.
   --tel tel.json             json file containing telescope related parameters.
   --ace ace.json             json file containing ace parameters.
   --ctl ctl.json             json file containing control parameters.
   --od outdir                name of the directory to put the outputs.
   --overwrite                if set, overwrite option activated.

"""

from docopt import docopt
import os
import json
import numpy as np
import astropy.io.ascii as asc
import astropy.io.fits as pf
from jis.photonsim.extract_json import mkDet, mkControlParams, mkTel
from jis.photonsim.wfe import wfe_model_z, calc_wfe
from jis.photonsim.response import calc_response
from jis.photonsim.psf import calc_psf
from jis.photonsim.ace import calc_ace


# Constants ########################################################
Rv = 3.1
JH = 2.0
alp = 0.75
spixdim  = [32, 32] # subpixel dimension in a pixel (setting for intrapix pattern).


# Command line interface
if __name__ == '__main__':
    args = docopt(__doc__)

    # Getting the parameters from command line. ####################
    dirname_params = ''
    if args['--pd']:
        dirname_params = args['--pd']

    filename_starplate = dirname_params + "/" + args['--starplate']
    filename_detjson   = dirname_params + "/" + args['--det']
    filename_teljson   = dirname_params + "/" + args['--tel']
    filename_acejson   = dirname_params + "/" + args['--ace']
    filename_ctljson   = dirname_params + "/" + args['--ctl']

    dirname_output = '.'
    if args['--od']:
        dirname_output = args['--od']

    overwrite = False
    if args['--overwrite']:
        overwrite = True


    # Setting output filenames. ####################################
    filename_interpix = dirname_output + "/" + "interpix.fits"
    filename_intrapix = dirname_output + "/" + "intrapix.fits"
    filename_wfejson  = dirname_output + "/" + "wfe.json"
    filename_wfe      = dirname_output + "/" + "wfe.fits"
    filename_aperture = dirname_output + "/" + "aperture.fits"
    filename_psf      = dirname_output + "/" + "psf.fits"
    filename_acex     = dirname_output + "/" + "aceX.fits"
    filename_acey     = dirname_output + "/" + "aceY.fits"

    filenames_output = [filename_interpix, filename_intrapix,\
                        filename_wfejson, filename_wfe,\
                        filename_aperture, filename_psf,\
                        filename_acex, filename_acey]


    # Checking the output directory. ###############################
    if not os.path.exists(dirname_output):
        os.makedirs(dirname_output)
    else:
        if overwrite is not True:
            for filename in filenames_output:
                if os.path.exists(filename):
                    print("\"{}\" exists.".format(filename))
                    print("Please set --overwrite option to overwrite it.")
                    exit()


    # Loading parameters. ##########################################
    table_starplate = asc.read(filename_starplate)
    detector        = mkDet(filename_detjson, spixdim=spixdim)
    control_params  = mkControlParams(filename_ctljson)
    telescope       = mkTel(filename_teljson)
    with open(filename_acejson, "r") as f:
        ace_params = json.load(f)
    f.close()
        

    # Selecting the data for the first plate. ######################
    pos = np.where(table_starplate['plate_id']==0)
    table_starplate = table_starplate[pos]


    # Making random wfe. ###########################################
    wp  = control_params.wfe_control
    wfe_amplitudes = wfe_model_z(np.random, wp['zernike_nmax'], wp['reference_wl'],\
                                 wp['zernike_odd'], wp['zernike_even'])

    # Saving amplitude data...
    with open(filename_wfejson, mode='w') as f:
        json.dump(wfe_amplitudes, f, indent=2)
    f.close()

    # Making wfe map...
    wfe = calc_wfe(telescope.epd, filename_wfejson)


    # Making PSFs ##################################################
    opteff = telescope.opt_efficiency 
    qe = detector.qe

    ## Currently, only one case of (Rv, JH).
    total_e_rate, wl_e_rate, e_rate =\
        calc_response(Rv, JH, alp,\
                      len(opteff['wl']), opteff['wl'], opteff['val'],\
                      np.min(opteff['wl']), np.max(opteff['wl']),\
                      qe['wl'], qe['val'])
    # total_e_rate in e/s/m^2; wl_e_rate in um; e_rate in e/s/m^2/um.
    # these values are for an object with an apparent Hw mag of 0 mag.

    ## Currently, only one PSF.
    print("Calculating PSF...")
    psf = calc_psf(wfe, wfe.shape[0],\
                   len(wl_e_rate), wl_e_rate, e_rate, total_e_rate,\
                   telescope.total_area, telescope.aperture,\
                   control_params.M_parameter, telescope.aperture.shape[0])
    # psf is that of an object which has the JH color of the set value and Hw=0.
    # The unit is e/sec/pix.

    # TK #######################################################################
    # Why do we need Rv and JH color excess?
    # I think we need apparent Hw mag and apparent J-H color instead of those.
    #
    # For considering various color objects, it might be good to calculate PSFs
    # with some J-H colors and use them with interpolating.
    # PSF calculation takes a long time.
    ############################################################################.


    # Ace simulation. ##############################################   
    ace_cp = control_params.ace_control
    print("Making ACE (X)...")
    acex, psdx = calc_ace(np.random, ace_cp['nace'], ace_cp['tace'], ace_params)

    print("Making ACE (Y)...")
    acey, psdy = calc_ace(np.random, ace_cp['nace'], ace_cp['tace'], ace_params)
 

    # Saving the outputs. ##########################################
    pf.writeto(filename_interpix, detector.interpix, overwrite=overwrite)
    pf.writeto(filename_intrapix, detector.intrapix, overwrite=overwrite)
    pf.writeto(filename_psf, psf, overwrite=overwrite)
    
    hdu = pf.PrimaryHDU(wfe)
    hdu.header["WFE-FILE"] = filename_wfejson
    hdu.header["WFE-EPD"]  = telescope.epd
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(filename_wfe, overwrite=overwrite)

    hdu = pf.PrimaryHDU(telescope.aperture)
    hdu.header["APTFILE"] = filename_teljson
    hdu.header["EPD"]     = telescope.epd
    hdu.header["COBS"]    = telescope.cobs
    hdu.header["STYPE"]   = telescope.spider_type
    hdu.header["STEL"]    = telescope.total_area * 1.e-6 # total area in m^2
    hdu.list = pf.HDUList([hdu])
    hdulist.writeto(filename_aperture, overwrite=overwrite)

    hdu = pf.PrimaryHDU(acex)
    hdu.header["ACE-FILE"] = filename_acejson
    hdu.header["ACE-TOTT"] = control_params.ace_control['tace']
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(filename_acex, overwrite=overwrite)

    hdu = pf.PrimaryHDU(acey)
    hdu.header["ACE-FILE"] = filename_acejson
    hdu.header["ACE-TOTT"] = control_params.ace_control['tace']
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(filename_acey, overwrite=overwrite)
