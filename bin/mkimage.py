#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""
  Make an image

  usage:
    mkimage.py [-h|--help] [--pd paramdir] --starplate star_plate.csv --det det.json 

 options:
   --help                     show this help message and exit
   --pd paramdir              name of the directory containing parameter files.
   --starplate star_plate.csv csv file containing star info (plate_id, star_id, xpix, ypix, l, b)
   --det det.json             json file containing detector related parameters.

"""

from docopt import docopt
import numpy as np
import astropy.io.ascii as asc
from jis.photonsim.extract_json import mkDet


# Command line interface
if __name__ == '__main__':
    args = docopt(__doc__)

    # Get parameters from command line
    dirname_params = ''
    if args['--pd']:
        dirname_params = args['--pd']
    filename_starplate = dirname_params + args['--starplate']
    filename_detjson   = dirname_params + args['--det']

    # Loading parameters.
    table_starplate = asc.read(filename_starplate)
    detector        = mkDet(filename_detjson)

    # Selecting the data for the first plate.
    pos = np.where(table_starplate['plate_id']==0)
    table_starplate = table_starplate[pos]

