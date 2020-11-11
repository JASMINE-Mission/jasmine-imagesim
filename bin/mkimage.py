#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""
  Make an image

  usage:
    mkimage.py [-h|--help] [--pd paramdir] --starplate star_plate.csv --det det.json [--od outdir] [--overwrite] 

 options:
   --help                     show this help message and exit.
   --pd paramdir              name of the directory containing parameter files.
   --starplate star_plate.csv csv file containing star info (plate_id, star_id, xpix, ypix, l, b)
   --det det.json             json file containing detector related parameters.
   --od outdir                name of the directory to put the outputs.
   --overwrite                if set, overwrite option activated.

"""

from docopt import docopt
import os
import numpy as np
import astropy.io.ascii as asc
import astropy.io.fits as pf
from jis.photonsim.extract_json import mkDet


# Command line interface
if __name__ == '__main__':
    args = docopt(__doc__)

    # Getting the parameters from command line. ####################
    dirname_params = ''
    if args['--pd']:
        dirname_params = args['--pd']

    filename_starplate = dirname_params + "/" + args['--starplate']
    filename_detjson   = dirname_params + "/" + args['--det']

    dirname_output = '.'
    if args['--od']:
        dirname_output = args['--od']

    overwrite = False
    if args['--overwrite']:
        overwrite = True


    # Setting output filenames. ####################################
    filename_interpix = dirname_output + "/" + "interpix.fits"


    # Checking the output directory. ###############################
    if not os.path.exists(dirname_output):
        os.makedirs(dirname_output)
    else:
        if overwrite is not True:
            for filename in [filename_interpix]:
                if os.path.exists(filename):
                    print("\"{}\" exists.".format(filename))
                    print("Please set --overwrite option to overwrite it.")
                    exit()


    # Loading parameters. ##########################################
    table_starplate = asc.read(filename_starplate)
    detector        = mkDet(filename_detjson)


    # Selecting the data for the first plate. ######################
    pos = np.where(table_starplate['plate_id']==0)
    table_starplate = table_starplate[pos]


    # Saving the outputs. ##########################################.
    pf.writeto(filename_interpix, detector.interpix, overwrite=overwrite)
