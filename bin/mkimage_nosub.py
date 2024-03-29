#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""Make an image.

.. code-block:: bash

  usage:
    mkimage_transit.py [-h|--help] [--pd paramdir] --starplate star_plate.csv [--var variability.json] --det det_planet.json --tel tel.json --ace ace.json --ctl ctl.json [--dft drift.json] --format format [--od outdir] [--overwrite]

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
from jis.photonsim.extract_json import Detector, ControlParams, Telescope, Variability, Drift
from jis.photonsim.wfe import wfe_model_z, calc_wfe, calc_dummy_wfe
from jis.photonsim.response import calc_response
from jis.photonsim.ace import calc_ace, calc_dummy_ace
from jis.photonsim.psf import calc_psf, calc_dummy_psf
from jis.pixsim import readflat as rf
from jis.pixsim import simpix_stable as sp
from jis.pixsim.integrate import integrate
from jis.pixsim.addnoise import addnoise
import matplotlib.pylab as plt

# Command line interface
if __name__ == '__main__':
    args = docopt(__doc__)

    # Getting the parameters from command line. ####################
    dirname_params = ''
    if args['--pd']:
        dirname_params = args['--pd']

    filename_starplate = os.path.join(dirname_params, args['--starplate'])

    if args['--var']:
        filename_varjson = os.path.join(dirname_params, args['--var'])
    filename_detjson = os.path.join(dirname_params, args['--det'])
    filename_teljson = os.path.join(dirname_params, args['--tel'])
    filename_acejson = os.path.join(dirname_params, args['--ace'])
    filename_ctljson = os.path.join(dirname_params, args['--ctl'])

    if args['--dft']:
        filename_dftjson = os.path.join(dirname_params, args['--dft'])

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

    # Loading parameters. ##########################################
    table_starplate = asc.read(filename_starplate)
    detector = Detector.from_json(filename_detjson)
    control_params = ControlParams.from_json(filename_ctljson)
    telescope = Telescope.from_json(filename_teljson)
    with open(filename_acejson, 'r') as f:
        ace_params = json.load(f)
    f.close()

    # det. pix. scale in arcsec/pix.
    detpix_scale = detector.pixsize*1.e-6/telescope.efl/1.e-3*180.*3600./np.pi

    # Setting output filenames. ####################################
    filename_interpix = os.path.join(dirname_output, 'interpix.fits')
    filename_intrapix = os.path.join(dirname_output, 'intrapix.fits')
    filename_wfejson = os.path.join(dirname_output, 'wfe.json')
    filename_wfe = os.path.join(dirname_output, 'wfe.fits')
    filename_aperture = os.path.join(dirname_output, 'aperture.fits')
    filename_psf = os.path.join(dirname_output, 'psf.fits')
    filename_acex = os.path.join(dirname_output, 'aceX.fits')
    filename_acey = os.path.join(dirname_output, 'aceY.fits')

    if output_format == 'platefits':
        filename_images = []
        for i in range(0, control_params.nplate):
            filename_images.append(os.path.join(
                dirname_output, 'image{:02d}.fits'.format(i)))
    elif output_format == 'fitscube':
        filename_images = [os.path.join(dirname_output, 'image.fits')]
    elif output_format == 'hdfcube':
        filename_images = [os.path.join(dirname_output, 'image.h5')]

    filenames_output = [filename_interpix, filename_intrapix,
                        filename_wfejson, filename_wfe,
                        filename_aperture, filename_psf,
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
                    print('Please set --overwrite option to overwrite it.')
                    exit()

    # Selecting the data for the first plate. ######################
    pos = np.where(table_starplate['plate index'] == 0)
    table_starplate = table_starplate[pos]

    # Making random wfe. ###########################################
    if control_params.effect.wfe is True:
        print('calculate WFE...')
        wp = control_params.wfe_control
        wfe_amplitudes = wfe_model_z(
            np.random, wp['zernike_nmax'], wp['reference_wl'],
            wp['zernike_odd'], wp['zernike_even'])

        # Saving amplitude data...
        with open(filename_wfejson, mode='w') as f:
            json.dump(wfe_amplitudes, f, indent=2)

        # Making wfe map...
        wfe = calc_wfe(telescope.epd, filename_wfejson)
    else:
        print('WFE simulation is skipped.')
        wfe = calc_dummy_wfe(telescope.epd, 'dummy.json')

    # Making PSFs ##################################################
    opteff = telescope.opt_efficiency
    qe = detector.qe

    # Currently, only one case of (Rv, JH).
    Rv = control_params.Rv
    JH = control_params.JH
    alpha = control_params.alpha
    total_e_rate, wl_e_rate, e_rate = \
        calc_response(Rv, JH, alpha, opteff.wavelength, opteff.efficiency,
                      np.min(opteff.wavelength), np.max(opteff.wavelength), qe.wl, qe.val)
    # total_e_rate in e/s/m^2; wl_e_rate in um; e_rate in e/s/m^2/um.
    # these values are for an object with an apparent Hw mag of 0 mag.

    if control_params.effect.psf is True:
        # Currently, only one PSF.
        print('Calculating PSF...')
        psf = calc_psf(wfe, wfe.shape[0],
                       wl_e_rate, e_rate, total_e_rate,
                       telescope.total_area, telescope.aperture,
                       control_params.M_parameter, telescope.aperture.shape[0],
                       control_params.fN_parameter)
        # psf is that of an object which has the JH color of the set value and Hw=0.
        # The unit is e/sec/pix.
    else:
        print('Realistic PSF simulation is skipped.')
        print('Generate fake PSF...')
        psf = calc_dummy_psf(wfe, wfe.shape[0],
                             wl_e_rate, e_rate, total_e_rate,
                             telescope.total_area, telescope.aperture,
                             control_params.M_parameter, telescope.aperture.shape[0])

    # Ace simulation. ##############################################
    nace = control_params.ace_control.get('nace')
    tace = control_params.ace_control.get('tace')
    if control_params.effect.ace is True:
        print('Making ACE (X)...')
        rg_acex = np.random.default_rng(
            control_params.ace_control.get('acex_seed'))
        acex, psdx = calc_ace(rg_acex, nace, tace, ace_params)
        # the standard deviation of acex is normalized to unity.

        print('Making ACE (Y)...')
        rg_acey = np.random.default_rng(
            control_params.ace_control.get('acey_seed'))
        acey, psdy = calc_ace(rg_acey, nace, tace, ace_params)
        # the standard deviation of acey is normalized to unity.
    else:
        print('ACE simulation is skipped.')
        print('Generate face ACE(X) and ACE(Y)...')
        acex = calc_dummy_ace(np.random, nace, tace, ace_params)
        acey = calc_dummy_ace(np.random, nace, tace, ace_params)

    tplate = control_params.tplate
    tscan = detector.readparams.t_scan
    dtace = control_params.ace_control['dtace']
    # Number of timesteps per a plate.
    Nts_per_plate = int((tplate+tscan)/dtace+0.5)

    # Drift
    if args['--dft']:
        dft = Drift.from_json(filename_dftjson)
        dft.compute_drift(dtace, nace)

    # Preparation for making image. ################################

    # Full data of the displacement in detpix.
    # (ace[x|y] scaled and converted to detpix)
    acex_std = control_params.ace_control.get('acex_std')
    acey_std = control_params.ace_control.get('acey_std')
    # Setting and plotting full trajectory.
    if args['--dft']:
        theta_full = np.array([acex*acex_std/detpix_scale+dft.drift_theta[0, :],
                               acey*acey_std/detpix_scale+dft.drift_theta[1, :]])
        plt.plot(acex*acex_std/detpix_scale +
                 dft.drift_theta[0, :], acey*acey_std/detpix_scale+dft.drift_theta[1, :], '.')
        plt.savefig('theta.png')
    else:
        theta_full = np.array(
            [acex*acex_std/detpix_scale, acey*acey_std/detpix_scale])

    Npixcube = int((np.max(np.abs(theta_full))+detector.nmargin)*2)
    pixdim = [Npixcube, Npixcube]  # adaptive pixel dimension in the aperture.

    # Variablity
    varsw = False
    if args['--var']:
        # load variability class
        variability = Variability.from_json(filename_varjson)
        # define time array in the unit of day
        tday = (tplate+tscan)*np.array(range(0, control_params.nplate))/3600/24
        for line in asc.read(filename_starplate):
            varsw, injlc, b = variability.read_var(tday, line['star index'])
            if varsw:
                plt.plot(tday, injlc)
                plt.savefig('variability_input'+'_' +
                            str(line['star index'])+'.png')
                plt.clf()

    if Nts_per_plate*control_params.nplate >= theta_full.shape[1]:
        print('Insufficient time length of ACE data.')
        print('Nts_per_plate*Nplate: {}'.format(Nts_per_plate*control_params.nplate))
        print('N_ace_data          : {}'.format(theta_full.shape[1]))
        sys.exit(-1)

    # psf center in the unit of fp-cell
    psfcenter = (np.array(np.shape(psf))-1.0)*0.5

    fp_cellsize_rad = (1./control_params.M_parameter)*1.e-3  # in rad/fp-cell.
    fp_scale = fp_cellsize_rad * 3600.*180./np.pi           # arcsec/fp-cell.
    psfscale = fp_scale/detpix_scale                        # det-pix/fp-cell.

    # Making image. ################################################
    uniform_flat_interpix = np.ones_like(detector.flat.interpix)
    uniform_flat_intrapix = np.ones_like(detector.flat.intrapix)

    # Making sky region.
    pixcube_global = np.zeros(
        shape=(detector.npix, detector.npix, control_params.nplate))
    pixcube_global += detector.idark * tplate
    pixcube_global, seed = addnoise(
        pixcube_global, np.sqrt(2.)*detector.readnoise)
    # in adu/pix/plate.
    pixcube_global = np.round(pixcube_global/detector.gain)

    pixcube_global_adu1 = np.zeros(
        shape=(detector.npix, detector.npix, control_params.nplate))
    pixcube_global_adu2 = np.zeros(
        shape=(detector.npix, detector.npix, control_params.nplate))

    # Making data around each star.
    for line in table_starplate:
        print('StarID: {}'.format(line['star index']))

        # Position setting.
        xc_global = line['x pixel'] - 1  # Stellar pos. in glob. coord (X).
        yc_global = line['y pixel'] - 1  # Stellar pos. in glob. coord (Y).
        # Origin pix position in global coord (x).
        x0_global = int(xc_global - Npixcube*0.5 + 0.5)
        # Origin pix position in global coord (y).
        y0_global = int(yc_global - Npixcube*0.5 + 0.5)
        xc_local = xc_global - x0_global  # Stellar position (local; x).
        yc_local = yc_global - y0_global  # Stellar position (local; y).
        mag = line['Hwmag']

        # Making local flat data.
        if control_params.effect.flat_interpix is True:
            interpix_local = rf.flat_interpix(
                detector.flat.interpix, x0_global, y0_global, pixdim, figsw=0)
        else:
            interpix_local = rf.flat_interpix(
                uniform_flat_interpix, x0_global, y0_global, pixdim, figsw=0)

        # Making a cube containing plate data for a local region (small frame for a local region).
        # Initialize (Axis order: X, Y, Z)
        pixcube = np.zeros((Npixcube, Npixcube, control_params.nplate))
        # Initialize (Axis order: X, Y, Z)
        pixcube_adu1 = np.zeros((Npixcube, Npixcube, control_params.nplate))
        # Initialize (Axis order: X, Y, Z)
        pixcube_adu2 = np.zeros((Npixcube, Npixcube, control_params.nplate))

        # Load variability
        if args['--var']:
            varsw, injlc, b = variability.read_var(tday, line['star index'])

        # Loop to take each plate.
        for iplate in tqdm.tqdm(range(0, control_params.nplate)):
            # picking temporary trajectory and local position update
            istart = iplate * Nts_per_plate
            iend = (iplate+1)*Nts_per_plate

            # Displacement from the initial position.
            theta = np.copy(theta_full[:, istart:iend])
            # Displacement in local coord.
            theta = theta + np.array([[xc_local, yc_local]]).T
            # 0.5-pix shift to treat the coodinate difference.
            theta = theta + np.array([[0.5, 0.5]]).T
            # Global coord.   : (0, 0) is the center of the bottom-left corner pixel.
            # Local coord.    : (0, 0) is the center of the bottom-left corner pixel..
            # Coord. in simpix: (0, 0) is the bottom-left corner of the bottom-left corner pixel.

            # Performing the PSF integration (simpix).
            #   Output: array of images in each time bin in the exposure.
            #   When the PSF is given in e/fp-cell/sec,
            #   simpix/(psfscale*psfscale) is in e/pix/(1./Nts_per_plate sec).
            if control_params.effect.flat_intrapix:
                pixar = sp.simpix(theta, interpix_local, detector.flat.intrapix,
                                  psfarr=psf, psfcenter=psfcenter, psfscale=psfscale)\
                    / (psfscale*psfscale)*dtace/(1./Nts_per_plate)
            else:
                pixar = sp.simpix(theta, interpix_local, uniform_flat_intrapix,
                                  psfarr=psf, psfcenter=psfcenter, psfscale=psfscale)\
                    / (psfscale*psfscale)*dtace/(1./Nts_per_plate)
            # pixar is in e/pix/dtace.

            # magnitude scaling.
            pixar = pixar * 10.**(mag/(-2.5))

            # variability
            """
            Curretly, the time resolution should be prepared in the unit of tplate + tscan. We do not support the finest time resolution yet (dtace).
            """
            if varsw:
                pixar = pixar*injlc[iplate]

            # Adding dark current (including stray light).
            dark = np.ones(shape=pixar.shape) * detector.idark * dtace
            pixar = pixar + dark

            # Integrating, adding noise, and quantization.
            adu2, adu1 = integrate(
                pixar, x0_global, y0_global, tplate, dtace, detector, raw=True)
            integrated = adu2 - adu1

            # integrated is in adu/pix/plate.
            pixcube[:, :, iplate] = integrated
            pixcube_adu1[:, :, iplate] = adu1
            pixcube_adu2[:, :, iplate] = adu2

            pixcube_global[x0_global:x0_global+Npixcube, y0_global:y0_global+Npixcube, iplate] =\
                pixcube[:, :, iplate]
            pixcube_global_adu1[x0_global:x0_global+Npixcube, y0_global:y0_global+Npixcube, iplate] =\
                pixcube_adu1[:, :, iplate]
            pixcube_global_adu2[x0_global:x0_global+Npixcube, y0_global:y0_global+Npixcube, iplate] =\
                pixcube_adu2[:, :, iplate]

    # Saving the outputs. ##########################################
    if control_params.effect.flat_interpix is True:
        pf.writeto(filename_interpix, detector.flat.interpix,
                   overwrite=overwrite)
    else:
        pf.writeto(filename_interpix, uniform_flat_interpix,
                   overwrite=overwrite)
    if control_params.effect.flat_interpix is True:
        pf.writeto(filename_intrapix, detector.flat.intrapix,
                   overwrite=overwrite)
    else:
        pf.writeto(filename_intrapix, uniform_flat_intrapix,
                   overwrite=overwrite)
    pf.writeto(filename_psf, psf, overwrite=overwrite)
    if output_format == 'hdfcube':
        with h5py.File(filename_images[0], 'w') as f:
            f.create_group('header')
            f.create_group('data')
            f.create_dataset('header/tplate', data=tplate)
            f.create_dataset('header/unit', data='e-/pix/plate')
            f.create_dataset('data/pixcube', data=pixcube_global)
            f.create_dataset('data/pixcube_adu1', data=pixcube_global_adu1)
            f.create_dataset('data/pixcube_adu2', data=pixcube_global_adu2)

    else:
        pixcube_global = np.swapaxes(pixcube_global, 0, 2)
        if output_format == 'platefits':
            for i in range(0, control_params.nplate):
                pf.writeto(filename_images[i], pixcube_global[i].astype(
                    'int32'), overwrite=overwrite)
        elif output_format == 'fitscube':
            pf.writeto(filename_images[0], pixcube_global.astype(
                'int32'), overwrite=overwrite)

    hdu = pf.PrimaryHDU(wfe)
    hdu.header['WFE-FILE'] = filename_wfejson
    hdu.header['WFE-EPD'] = telescope.epd
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(filename_wfe, overwrite=overwrite)

    hdu = pf.PrimaryHDU(telescope.aperture)
    hdu.header['APTFILE'] = filename_teljson
    hdu.header['EPD'] = telescope.epd
    hdu.header['COBS'] = telescope.cobs
    hdu.header['STYPE'] = telescope.spider.type
    hdu.header['STEL'] = telescope.total_area  # total area in m^2
    hdu.list = pf.HDUList([hdu])
    hdulist.writeto(filename_aperture, overwrite=overwrite)

    hdu = pf.PrimaryHDU(acex)
    hdu.header['ACE-FILE'] = filename_acejson
    hdu.header['ACE-TOTT'] = control_params.ace_control['tace']
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(filename_acex, overwrite=overwrite)

    hdu = pf.PrimaryHDU(acey)
    hdu.header['ACE-FILE'] = filename_acejson
    hdu.header['ACE-TOTT'] = control_params.ace_control['tace']
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(filename_acey, overwrite=overwrite)
