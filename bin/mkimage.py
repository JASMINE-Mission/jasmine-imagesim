#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""
  Make an image

  usage:
    mkimage.py [-h|--help] [--pd paramdir] --starplate star_plate.csv --det det.json --tel tel.json --ctl ctl.json [--od outdir] [--overwrite] 

 options:
   --help                     show this help message and exit.
   --pd paramdir              name of the directory containing parameter files.
   --starplate star_plate.csv csv file containing star info (plate_id, star_id, xpix, ypix, l, b)
   --det det.json             json file containing detector related parameters.
   --tel tel.json             json file containing telescope related parameters.
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

    filenames_output = [filename_interpix, filename_intrapix,\
                        filename_wfejson, filename_wfe]


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
    detector        = mkDet(filename_detjson)
    control_params  = mkControlParams(filename_ctljson)
    telescope       = mkTel(filename_teljson)


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

    # Making wfe map...
    wfe = calc_wfe(telescope.epd, filename_wfejson)


    # Saving the outputs. ##########################################.
    pf.writeto(filename_interpix, detector.interpix, overwrite=overwrite)
    pf.writeto(filename_intrapix, detector.intrapix, overwrite=overwrite)
    
    hdu = pf.PrimaryHDU(wfe)
    hdu.header["WFE-FILE"] = filename_wfejson
    hdu.header["WFE-EPD"]  = telescope.epd
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(filename_wfe, overwrite=overwrite)