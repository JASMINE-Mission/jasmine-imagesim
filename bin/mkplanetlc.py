#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""Make a planet signal 

.. code-block:: bash

  usage:
    mkplanetlc.py [-h|--help] [-m|--plot] -s tframe -f nframe -p planet.json -o lc.fits

  options:
    -h --help       show this help message and exit
    -m --plot       plot light curve
    -s tframe       exposure [sec] of a frame
    -f nframe       number of the frames
    -p planet.json  planet parameters
    -o lc.fits      output lightcurve fits
""" 
from docopt import docopt             # command line interface
import numpy as np
from jis.pixsim import transitmodel
import astropy.io.fits as fits

# Command line interface
if __name__ == '__main__':
    args = docopt(__doc__)
    # Get parameters from command line
    tframe=float(args['-s'])
    nframe=int(args['-f'])
    t=np.array(range(0,nframe))*tframe #total time sec
    t_day=t/3600/24
    
    #planet model
    Rpin=1.0
    injlc, b=transitmodel.gentransit_json(t_day,args['-p'])

    #plotting
    if args['--plot']:
        import matplotlib.pyplot as plt
        plt.plot(t,injlc)
        plt.xlabel("Time sec")
        plt.savefig("lcin.png")
        plt.show()

    # Save Light Curve
    data=np.array([t,injlc])
    hdu = fits.PrimaryHDU(data)
    hdu.header["EXPFRAME"] = args['-s']
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(args['-o'],overwrite=True)
