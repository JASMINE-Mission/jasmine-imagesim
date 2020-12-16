#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""Make an image

  usage:
    mkimage.py [-h|--help] [--pd paramdir] --starplate star_plate.csv [--var variability.json] --det det.json --tel tel.json --ace ace.json --ctl ctl.json --format format [--od outdir] [--overwrite] 

  options:
   --help                     show this help message and exit.
   --pd paramdir              name of the directory containing parameter files.
   --starplate star_plate.csv csv file containing star info (plate index, star index, x pixel, y pixel, lambda, beta)
   --var variability.json     json file for stellar variability/transit (optional). The input variability will be shown in variability_input().png 
   --det det.json             json file containing detector related parameters.
   --tel tel.json             json file containing telescope related parameters.
   --ace ace.json             json file containing ace parameters.
   --ctl ctl.json             json file containing control parameters.
   --format format            format of the output file (platefits, fitscube, hdfcube).
   --od outdir                name of the directory to put the outputs.
   --overwrite                if set, overwrite option activated.

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
from jis.photonsim.extract_json import mkDet, mkControlParams, mkTel, mkVar
from jis.photonsim.wfe import wfe_model_z, calc_wfe
from jis.photonsim.response import calc_response
from jis.photonsim.psf import calc_psf
from jis.photonsim.ace import calc_ace
from jis.pixsim import readflat as rf
from jis.pixsim import simpix_stable as sp
from jis.pixsim.integrate import integrate
from jis.pixsim.addnoise import addnoise
import matplotlib.pylab as plt


# Constants ########################################################
Rv  = 3.1
JH  = 2.0
alp = 0.75
spixdim  = [32, 32]  # subpixel dimension in a pixel (setting for intrapix pattern).
acex_std = 0.276     # std of ace_x (arcsec).
acey_std = 0.276     # std of ace_y (arcsec). 
Nmargin = 10         # Margin for simpix calc.
Nplate  = 200         # Number of plates in a small frame.
tplate  = 12.5       # Exposure time of a plate (sec).
mag     = 12.5       # Stellar Hw band mag.


# Command line interface
if __name__ == '__main__':
    args = docopt(__doc__)

    # Getting the parameters from command line. ####################
    dirname_params = ''
    if args['--pd']:
        dirname_params = args['--pd']

    filename_starplate = os.path.join(dirname_params, args['--starplate'])
    filename_varjson   = os.path.join(dirname_params, args['--var'])
    filename_detjson   = os.path.join(dirname_params, args['--det'])
    filename_teljson   = os.path.join(dirname_params, args['--tel'])
    filename_acejson   = os.path.join(dirname_params, args['--ace'])
    filename_ctljson   = os.path.join(dirname_params, args['--ctl'])

    varsw=False
    if args['--var']:
        #load variability class
        variability=mkVar(filename_varjson)
        #define time array in the unit of day
        tday=tplate*np.array(range(0,Nplate))/3600/24
        for line in asc.read(filename_starplate):
            varsw, injlc, b=variability.read_var(tday,line['plate index'],line['star index'])
            if varsw:
                plt.plot(tday,injlc)
                plt.savefig("variability_input"+str(line['plate index'])+"_"+str(line['star index'])+".png")
                plt.clf()
                
    output_format = args['--format']
    if output_format not in ['platefits', 'fitscube', 'hdfcube']:
        print("format must be 'platefits', 'fitscube' or 'hdfcube'.")
        exit(-1)

    dirname_output = '.'
    if args['--od']:
        dirname_output = args['--od']

    overwrite = False
    if args['--overwrite']:
        overwrite = True


    # Setting output filenames. ####################################
    filename_interpix = os.path.join(dirname_output, "interpix.fits")
    filename_intrapix = os.path.join(dirname_output, "intrapix.fits")
    filename_wfejson  = os.path.join(dirname_output, "wfe.json")
    filename_wfe      = os.path.join(dirname_output, "wfe.fits")
    filename_aperture = os.path.join(dirname_output, "aperture.fits")
    filename_psf      = os.path.join(dirname_output, "psf.fits")
    filename_acex     = os.path.join(dirname_output, "aceX.fits")
    filename_acey     = os.path.join(dirname_output, "aceY.fits")

    if output_format == 'platefits':
        filename_images = []
        for i in range(0, Nplate):
            filename_images.append(os.path.join(dirname_output, "image{:02d}.fits".format(i)))
    elif output_format == 'fitscube':
        filename_images = [os.path.join(dirname_output, "image.fits")]
    elif output_format == 'hdfcube':
        filename_images = [os.path.join(dirname_output, "image.h5")]

    filenames_output = [filename_interpix, filename_intrapix,\
                        filename_wfejson, filename_wfe,\
                        filename_aperture, filename_psf,\
                        filename_acex, filename_acey]
    filenames_output = filenames_output + filename_images


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

    detpix_scale = detector.pixsize*1.e-6/telescope.efl/1.e-3*180.*3600./np.pi # det. pix. scale in arcsec/pix.
        

    # Selecting the data for the first plate. ######################
    pos = np.where(table_starplate['plate index']==0)
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
    # acex is normalized by the std.

    print("Making ACE (Y)...")
    acey, psdy = calc_ace(np.random, ace_cp['nace'], ace_cp['tace'], ace_params)
    # acey is normalized by the std.


    # Preparation for making image. ################################

    ## Full data of the displacement in detpix.
    ## (ace[x|y] scaled and converted to detpix)
    theta_full = np.array([acex*acex_std/detpix_scale, acey*acey_std/detpix_scale])

    Npixcube = int((np.max(np.abs(theta_full))+Nmargin)*2)
    pixdim   = [Npixcube, Npixcube] # adaptive pixel dimension in the aperture.

    tscan = detector.t_overhead +\
            detector.tsmpl*(detector.npix_pre+detector.ncol_ch+detector.npix_post)*detector.nrow_ch
    dtace = control_params.ace_control['dtace']
    Nts_per_plate = int((tplate+tscan)/dtace+0.5) # Number of timesteps per a plate.

    if Nts_per_plate*Nplate >= theta_full.shape[1]:
        print("Insufficient time length of ACE data.")
        print("Nts_per_plate*Nplate: {}".format(Nts_per_plate*Nplate))
        print("N_ace_data          : {}".format(theta_full.shape[1]))
        sys.exit(-1)

    psfcenter = (np.array(np.shape(psf))-1.0)*0.5 #psf center in the unit of fp-cell

    fp_cellsize_rad = (1./control_params.M_parameter)*1.e-3 # in rad/fp-cell.
    fp_scale = fp_cellsize_rad * 3600.*180./np.pi           # arcsec/fp-cell.
    psfscale = fp_scale/detpix_scale                        # det-pix/fp-cell.


    # Making image. ################################################

    ## Making sky region.
    pixcube_global = np.zeros(shape=(detector.npix, detector.npix, Nplate))
    pixcube_global += detector.idark * tplate
    pixcube_global, seed = addnoise(pixcube_global, np.sqrt(2.)*detector.readnoise)
    pixcube_global = np.round(pixcube_global/detector.gain) # in adu/pix/plate.

    ## Making data around each star.
    for line in table_starplate:
        print("StarID: {}".format(line['star index']))

        # Position setting.
        xc_global = line['x pixel'] - 1 # Stellar pos. in glob. coord (X).
        yc_global = line['y pixel'] - 1 # Stellar pos. in glob. coord (Y).
        x0_global = int(xc_global - Npixcube*0.5 + 0.5) # Origin pix position in global coord (x).
        y0_global = int(yc_global - Npixcube*0.5 + 0.5) # Origin pix position in global coord (y).
        xc_local  = xc_global - x0_global  # Stellar position (local; x).
        yc_local  = yc_global - y0_global  # Stellar position (local; y).

        # Making local flat data.
        interpix_local = rf.flat_interpix(detector.interpix, x0_global, y0_global, pixdim, figsw=0)

        # Making a cube containing plate data for a local region (small frame for a local region).
        pixcube = np.zeros((Npixcube, Npixcube, Nplate))   # Initialize (Axis order: X, Y, Z)

        # Load variability
        if args["--var"]:
            varsw, injlc, b=variability.read_var(tday,line['plate index'],line['star index'])
        
        for iplate in tqdm.tqdm(range(0, Nplate)):         # Loop to take each plate.
            # picking temporary trajectory and local position update
            istart = iplate    *Nts_per_plate
            iend   = (iplate+1)*Nts_per_plate

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
            pixar = sp.simpix(theta, interpix_local, detector.intrapix,\
                              psfarr=psf, psfcenter=psfcenter, psfscale=psfscale)\
                              /(psfscale*psfscale)*dtace/(1./Nts_per_plate)
            # pixar is in e/pix/dtace.

            # magnitude scaling.
            pixar = pixar * 10.**(mag/(-2.5))

            # variability
            if varsw:
                pixar=pixar*injlc[iplate]

            # Adding dark current (including stray light).
            dark  = np.ones(shape=pixar.shape) * detector.idark * dtace
            pixar = pixar + dark

            # Integrating, adding noise, and quantization.
            integrated = integrate(pixar, x0_global, y0_global, tplate, dtace, detector)
            # integrated is in adu/pix/plate.
            pixcube[:,:,iplate] = integrated
            
            pixcube_global[x0_global:x0_global+Npixcube, y0_global:y0_global+Npixcube, iplate] =\
                pixcube[:,:,iplate]
 

    # Saving the outputs. ##########################################
    pf.writeto(filename_interpix, detector.interpix, overwrite=overwrite)
    pf.writeto(filename_intrapix, detector.intrapix, overwrite=overwrite)
    pf.writeto(filename_psf, psf, overwrite=overwrite)
    if output_format == 'hdfcube':
        with h5py.File(filename_images[0],"w") as f:
            f.create_group("header")
            f.create_group("data")
            f.create_dataset("header/tplate", data=tplate)
            f.create_dataset("header/unit", data="e-/pix/plate")
            f.create_dataset("data/pixcube", data=pixcube_global)
    else:
        pixcube_global = np.swapaxes(pixcube_global, 0, 2)
        if output_format == 'platefits':
            for i in range(0, Nplate):
                pf.writeto(filename_images[i], pixcube_global[i].astype('int32'), overwrite=overwrite)
        elif output_format == 'fitscube':
            pf.writeto(filename_images[0], pixcube_global.astype('int32'), overwrite=overwrite)
    
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
    hdu.header["STEL"]    = telescope.total_area  # total area in m^2
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
    
