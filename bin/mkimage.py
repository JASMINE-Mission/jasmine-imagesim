#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""Make an image

.. code-block:: bash

  usage:
    mkimage.py [-h|--help] [--pd paramdir] --starplate star_plate.csv [--var variability.json] --det det.json --tel tel.json --ace ace.json --ctl ctl.json [--dft drift.json] --format format [--od outdir] [--overwrite]

  options:
   -h --help                   show this help message and exit.
   --pd paramdir               name of the directory containing parameter files.
   --starplate star_plate.csv  csv file containing star info (plate index, star index, x pixel, y pixel, lambda, beta, Hwmag)
   --var variability.json      json file for stellar variability/transit (optional). The input variability will be shown in variability_input().png
   --det det.json              json file containing detector related parameters.
   --tel tel.json              json file containing telescope related parameters.
   --ace ace.json              json file containing ace parameters.
   --ctl ctl.json              json file containing control parameters.
   --dft drift.json            json file containing drift parameterers.
   --format format             format of the output file (platefits, fitscube, hdfcube).
   --od outdir                 name of the directory to put the outputs.
   --overwrite                 if set, overwrite option activated.

"""

from docopt import docopt
import os
import sys
import json
import tqdm
import h5py
import numpy as np
import astropy.io.ascii as asc
import astropy.io.fits as pf
from jis.binutils.setfiles import set_filenames_from_args, set_filenames_output, set_output_from_args, check_output_directory
from jis.binutils.setcontrol import load_parameters
from jis.binutils.save import save_outputs
from jis.binutils.run import run_wfe, run_psf, run_ace, init_pix

from jis.photonsim.extract_json import Detector, ControlParams, Telescope, Variability, Drift
from jis.pixsim import readflat as rf
from jis.pixsim import simpix_stable as sp
from jis.pixsim.integrate import integrate
from jis.pixsim.addnoise import addnoise
from scipy import ndimage
import matplotlib.pylab as plt

def get_detpix_scale(telescope,detector):
    """ det. pix. scale in arcsec/pix.

    Returns:
        detector pixel scale

    """
    return detector.pixsize*1.e-6/telescope.efl/1.e-3*180.*3600./np.pi 

def include_variability(control_params,detector):
    """ Variablity

    Args:
        

    """
    from jis.photonsim.extract_json import Variability

    varsw=False
    #load variability class
    variability=Variability.from_json(filenames["varjson"])
    #define time array in the unit of day
    tday=(control_params.tplate+detector.readparams.t_scan)*np.array(range(0,control_params.nplate))/3600/24
    for line in asc.read(filenames["starplate"]):
        varsw, injlc, b=variability.read_var(tday,line['star index'])
        if varsw:
            plt.plot(tday,injlc)
            plt.savefig("variability_input"+"_"+str(line['star index'])+".png")
            plt.clf()
    return variability, tday


# Command line interface
if __name__ == '__main__':
    
    args = docopt(__doc__)        
    filenames, dirname_output=set_filenames_from_args(args)
    table_starplate, detector, control_params, telescope, ace_params = load_parameters(filenames)
    filenames, output_format, overwrite=set_filenames_output(args,filenames,control_params,dirname_output)

    # Selecting the data for the first plate. ######################
    pos = np.where(table_starplate['plate index']==0)
    table_starplate = table_starplate[pos]
    
    wfe=run_wfe(control_params, telescope)
    psf=run_psf(control_params, telescope, detector, wfe)
    acex, acey, Nts_per_plate=run_ace(control_params, detector, ace_params)

    detpix_scale = get_detpix_scale(telescope,detector)
    theta_full, pixdim, Npixcube = init_pix(control_params,detector,acex,acey, detpix_scale,args["--dft"])
    
    if args['--var']:
        variability, tday = include_variability(control_params,detector)
        
    if Nts_per_plate*control_params.nplate >= theta_full.shape[1]:
        print("Insufficient time length of ACE data.")
        print("Nts_per_plate*Nplate: {}".format(Nts_per_plate*control_params.nplate))
        print("N_ace_data          : {}".format(theta_full.shape[1]))
        sys.exit(-1)

    fp_cellsize_rad = (1./control_params.M_parameter)*1.e-3 # in rad/fp-cell.
    fp_scale = fp_cellsize_rad * 3600.*180./np.pi           # arcsec/fp-cell.
    psfscale = fp_scale/detpix_scale                        # det-pix/fp-cell.

    ## In the gauss-ace mode, apply gauss filter to psf, here.
    if control_params.effect.ace == "gauss":
        if acex_std != acey_std:
            print("In the current gauss-ace mode, acex_std must be equal to acey_std. Sorry!")
            exit(-1)
        else:
            psf = ndimage.gaussian_filter(psf, sigma=acex_std/fp_scale)


    # Making image. ################################################
    uniform_flat_interpix = np.ones_like(detector.flat.interpix)
    uniform_flat_intrapix = np.ones_like(detector.flat.intrapix)

    ## Making sky region.
    pixcube_global = np.zeros(shape=(detector.npix, detector.npix, control_params.nplate))
    pixcube_global += detector.idark * control_params.tplate
    pixcube_global, seed = addnoise(pixcube_global, np.sqrt(2.)*detector.readnoise)
    pixcube_global = np.round(pixcube_global/detector.gain) # in adu/pix/plate.

    ## Making data around each star.
    for i_star in range(np.size(table_starplate)):
        line = table_starplate[i_star]
        print("StarID: {}".format(line['star index']))

        # Position setting.
        xc_global = line['x pixel'] - 1 # Stellar pos. in glob. coord (X).
        yc_global = line['y pixel'] - 1 # Stellar pos. in glob. coord (Y).
        x0_global = int(xc_global - Npixcube*0.5 + 0.5) # Origin pix position in global coord (x).
        y0_global = int(yc_global - Npixcube*0.5 + 0.5) # Origin pix position in global coord (y).
        xc_local  = xc_global - x0_global  # Stellar position (local; x).
        yc_local  = yc_global - y0_global  # Stellar position (local; y).
        mag = line['Hwmag']

        # Making local flat data.
        if control_params.effect.flat_interpix is True:
            interpix_local = rf.flat_interpix(
                detector.flat.interpix, x0_global, y0_global, pixdim, figsw=0)
        else:
            interpix_local = rf.flat_interpix(
                uniform_flat_interpix, x0_global, y0_global, pixdim, figsw=0)

        # Making a cube containing plate data for a local region (small frame for a local region).
        pixcube = np.zeros((Npixcube, Npixcube, control_params.nplate))   # Initialize (Axis order: X, Y, Z)

        # Load variability
        if args["--var"]:
            varsw, injlc, b=variability.read_var(tday,line['star index'])

        for iplate in tqdm.tqdm(range(0, control_params.nplate)):         # Loop to take each plate.
            # picking temporary trajectory and local position update
            istart = iplate    *Nts_per_plate
            iend   = (iplate+1)*Nts_per_plate
            # In no-ace mode, we make a single image with simpix to reduce the calculation time.
            # Below is a trick for that. After executing simpix, we will copy it to make Nts_per_plate shots.
            if control_params.effect.ace != "real":
                iend=istart+1
            
            theta = np.copy(theta_full[:,istart:iend])         # Displacement from the initial position.
            theta = theta + np.array([[xc_local, yc_local]]).T # Displacement in local coord.
            theta = theta + np.array([[0.5, 0.5]]).T           # 0.5-pix shift to treat the coodinate difference.
            # Global coord.   : (0, 0) is the center of the bottom-left corner pixel.
            # Local coord.    : (0, 0) is the center of the bottom-left corner pixel..
            # Coord. in simpix: (0, 0) is the bottom-left corner of the bottom-left corner pixel.

            # Performing the PSF integration (simpix).
            #   Output: array of images in each time bin in the exposure.
            #   When the PSF is given in e/fp-cell/sec,
            #   simpix/(psfscale*psfscale) is in e/pix/(1./Nts_per_plate sec).

            if control_params.effect.flat_intrapix:
                flat_intrapix = detector.flat.intrapix
            else:
                flat_intrapix = uniform_flat_intrapix

            if control_params.effect.wfe != 'fringe37':
                psfcenter = (np.array(np.shape(psf))-1.0)*0.5 #psf center in the unit of fp-cell
                pixar = sp.simpix(theta, interpix_local, flat_intrapix,\
                                  psfarr=psf, psfcenter=psfcenter, psfscale=psfscale)\
                                  /(psfscale*psfscale)*control_params.ace_control['dtace']/(1./Nts_per_plate)
            else:
                psfcenter = (np.array(np.shape(psf)[1:])-1.0)*0.5 #psf center in the unit of fp-cell
                pixar = sp.simpix(theta, interpix_local, flat_intrapix,\
                                  psfarr=psf[i_star], psfcenter=psfcenter, psfscale=psfscale)\
                                  /(psfscale*psfscale)*control_params.ace_control['dtace']/(1./Nts_per_plate)
            # pixar is in e/pix/control_params.ace_control['dtace'].

            # In none/gauss mode, we copy the single-shot image to make the full-movie cube.
            if control_params.effect.ace != "real":
                upixar=pixar[:,:,0]
                nxt,nyt=np.shape(upixar)
                pixar=upixar[:,:,np.newaxis]+np.zeros((nxt,nyt,Nts_per_plate))
                pixar=pixar/Nts_per_plate
                # In the above process to make pixar, Nts_per_plate is multiplied
                # to the result of simpix to make the units of pixar to be e/pix/control_params.ace_control['dtace'].
                # But, in none/gauss mode, the scaling is not correct for simulating a single-shot image.
                # Therefore, we divide pixar by Nts_per_plate for correction.


            # magnitude scaling.
            pixar = pixar * 10.**(mag/(-2.5))

            # variability
            """
            Curretly, the time resolution should be prepared in the unit of control_params.tplate + detector.readparams.t_scan. We do not support the finest time resolution yet (control_params.ace_control['dtace']).
            """
            if varsw:
                pixar=pixar*injlc[iplate]

            # Adding dark current (including stray light).
            dark  = np.ones(shape=pixar.shape) * detector.idark * control_params.ace_control['dtace']
            pixar = pixar + dark

            # Integrating, adding noise, and quantization.
            integrated = integrate(pixar, x0_global, y0_global, control_params.tplate, control_params.ace_control['dtace'], detector)
            # integrated is in adu/pix/plate.
            pixcube[:,:,iplate] = integrated

            pixcube_global[x0_global:x0_global+Npixcube, y0_global:y0_global+Npixcube, iplate] =\
                pixcube[:,:,iplate]


    save_outputs(filenames, output_format, control_params, telescope, detector, wfe, psf, pixcube_global, control_params.tplate, uniform_flat_interpix, uniform_flat_intrapix, acex, acey, overwrite)
