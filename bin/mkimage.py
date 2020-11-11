#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""
  Make an image

  usage:
    mkimage.py [-h|--help] --starplate star_plate.csv 

 options:
   --help                     show this help message and exit
   --starplate star_plate.csv csv file containing star info (plate_id, star_id, xpix, ypix, l, b)
"""

from docopt import docopt
import numpy as np
import astropy.io.ascii as asc


# Command line interface
if __name__ == '__main__':
    args = docopt(__doc__)

    # Get parameters from command line
    filename_starplate = args['--starplate']

    # Loading parameters.
    table_starplate = asc.read(filename_starplate)

    # Selecting the data for the first plate.
    pos = np.where(table_starplate['plate_id']==0)
    table_starplate = table_starplate[pos]

    print(table_starplate)
