#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""Make a galactic center image.

"""

import tqdm
import numpy as np
from jis.binutils.save import save_outputs
from jis.binutils.runphotonsim import run_calc_wfe, run_calc_psf, run_calc_ace
from jis.binutils.runpixsim import init_pix, uniform_flat, init_images, set_positions, make_local_flat, index_control_trajectory, calc_theta, scaling_pixar, add_varability, add_dark_current, run_simpix
from jis.binutils.scales import get_pixelscales
from jis.binutils.check import check_ace_length
from jis.pixsim.integrate import integrate
import astropy.io.ascii as asc  
import matplotlib.pylab as plt
from jis.pixsim.wcs import set_wcs
from jis.photonsim.extract_json import Detector, ControlParams, Telescope, AcePsd
from jis.galcen.read_galcen_position import load_jscon_random_stars
from jis.galcen.read_galcen_position import random_stars_to_starplate
from jis.galcen.read_galcen_position import maximum_separation    

if __name__ == '__main__':
    """
    Notes:
        When performing the PSF integration (simpix),
        Output= array of images in each time bin in the exposure.
        When the PSF is given in e/fp-cell/sec,
        simpix/(psfscale*psfscale) is in e/pix/(1./Nts_per_plate sec).

    """
    import os
    dirname_params = "../params/galcenimg"
    filenames = {}
    #filenames['starplate'] = os.path.join(dirname_params, "gelcen_random_stars.json")
    filenames['starplate'] = os.path.join(dirname_params,
                                          "gelcen_random_stars.json")
    filenames['detjson'] = os.path.join(dirname_params, "det.json")
    filenames['teljson'] = os.path.join(dirname_params, "tel.json")
    filenames['acejson'] = os.path.join(dirname_params, "ace_001.json")
    filenames['ctljson'] = os.path.join(dirname_params, "ctl.json")

    detector = Detector.from_json(filenames['detjson'])
    control_params = ControlParams.from_json(filenames['ctljson'])
    telescope = Telescope.from_json(filenames['teljson'])
    ace_params = AcePsd.from_json(filenames['acejson']).parameters


    # Load Aizawa's random stars generated by jscon and convert it to the star plate form
    random_star_data = load_jscon_random_stars()
    maxsep = maximum_separation(random_star_data)
    print("maxsep=",maxsep,"deg")
    
    telescope_center_ra = np.median(random_star_data["ra"])
    telescope_center_dec = np.median(random_star_data["dec"])
    jasmine_wcs = set_wcs(telescope_center_ra, telescope_center_dec, detector,
                          telescope)
    coord = SkyCoord(ramin * u.degree,
                     decmin * u.degree)
    xpixel, ypixel = jasmine_wcs.world_to_pixel(coord)
    print(xpixel,ypixel,"-")

    table_starplate = random_stars_to_starplate(random_star_data, jasmine_wcs)
    print(np.min(table_starplate["x pixel"]))
    print(np.min(table_starplate["y pixel"]))
    print(np.max(table_starplate["x pixel"]))
    print(np.max(table_starplate["y pixel"]))
    
    import sys
    sys.exit()

    # Setting scaling constants. ###################################
    detpix_scale, fp_cellsize_rad, fp_scale, psfscale = get_pixelscales(
        control_params, telescope, detector)

    # Running calculations. ########################################
    wfe = run_calc_wfe(control_params, telescope, filenames)
    psf = run_calc_psf(control_params, telescope, detector, wfe)
    acex, acey, Nts_per_plate = run_calc_ace(control_params, detector,
                                             ace_params)
    # acex and acey are normalized with their stddevs.
    driftsw = False
    theta_full, pixdim, Npixcube = init_pix(filenames, control_params,
                                            detector, acex, acey, detpix_scale,
                                            driftsw)
    check_ace_length(Nts_per_plate, control_params, theta_full)

    uniform_flat_interpix, uniform_flat_intrapix = uniform_flat(detector)
    pixcube_global = init_images(control_params, detector)


    # Making data around each star.
    for i_star, line in enumerate(table_starplate):
        print(i_star,line)
        print(type(line['star index']))
        
        print('StarID: {}'.format(line['star index']))
        mag = line['Hwmag']
        xc_local, yc_local, x0_global, y0_global, xc_global, yc_global = set_positions(
            line, Npixcube)
        interpix_local = make_local_flat(control_params, detector, x0_global,
                                         y0_global, pixdim)

        # Making a cube containing plate data for a local region (small frame for a local region).
        # Initialize (Axis order: X, Y, Z)
        pixcube = np.zeros((Npixcube, Npixcube, control_params.nplate))


        # Loop to take each plate.
        for iplate in tqdm.tqdm(range(0, control_params.nplate)):
            # picking temporary trajectory and local position update
            istart, iend = index_control_trajectory(control_params, iplate,
                                                    Nts_per_plate)
            theta = calc_theta(theta_full, istart, iend, xc_local, yc_local)

            if control_params.effect.flat_intrapix:
                flat_intrapix = detector.flat.intrapix
            else:
                flat_intrapix = uniform_flat_intrapix

            if control_params.effect.wfe != 'fringe37':
                psfin = psf
                psfcenter = (np.array(np.shape(psfin)) - 1.0) * 0.5
            else:
                psfin = psf[i_star]
                psfcenter = (np.array(np.shape(psfin)[1:]) - 1.0) * 0.5

            pixar = run_simpix(control_params, theta, interpix_local,
                               flat_intrapix, psfin, psfcenter, psfscale,
                               Nts_per_plate)
            pixar = scaling_pixar(pixar, mag)


            pixar = add_dark_current(control_params, detector, pixar)
            integrated = integrate(pixar, x0_global, y0_global,
                                   control_params.tplate,
                                   control_params.ace_control['dtace'],
                                   detector)
            # integrated is in adu/pix/plate.
            pixcube[:, :, iplate] = integrated
            pixcube_global[x0_global:x0_global+Npixcube, y0_global:y0_global+Npixcube, iplate] =\
                pixcube[:, :, iplate]

#    save_outputs(filenames, output_format, control_params, telescope, detector,
#                 wfe, psf, pixcube_global, control_params.tplate,
#                 uniform_flat_interpix, uniform_flat_intrapix, acex, acey,
#                 overwrite)
