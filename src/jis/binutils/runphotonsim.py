import numpy as np
import json
from scipy import ndimage
from jis.photonsim.wfe import wfe_model_z, calc_wfe, calc_dummy_wfe, calc_wfe_fringe37
from jis.photonsim.psf import calc_psf, calc_gauss_psf
from jis.photonsim.response import calc_response
from jis.photonsim.ace import calc_ace, calc_dummy_ace


def run_calc_wfe(control_params, telescope, filenames):
    """Making wfe.

    Args:
        control_params: control parameters
        telescope: telescope object
    Returns:
        wavefront error
    """
    if control_params.effect.wfe == 'dummy':
        print('WFE simulation is skipped (making dummy).')
        wfe = calc_dummy_wfe(telescope.epd, 'dummy.json')
    elif control_params.effect.wfe == 'random':
        print('calculate random WFE...')
        wp = control_params.wfe_control
        wfe_amplitudes = wfe_model_z(
            np.random, wp['zernike_nmax'], wp['reference_wl'],
            wp['zernike_odd'], wp['zernike_even'])

        # Saving amplitude data...
        with open(filenames['wfejson'], mode='w') as f:
            json.dump(wfe_amplitudes, f, indent=2)

        # Making wfe map...
        wfe = calc_wfe(telescope.epd, filenames['wfejson'])
    elif control_params.effect.wfe == 'fringe37':
        print('calculate WFE with fringe37 params...')
        wp = control_params.wfe_control

        # Making position array (in deg).
        positions = np.array([table_starplate['x pixel']-1.+detector.offset_x_mm/detector.pixsize/1.e-3,
                              table_starplate['y pixel']-1.+detector.offset_y_mm/detector.pixsize/1.e-3]).T\
            * detpix_scale/3600.
        # detector.offset_[x|y]_mm is the position of (0, 0) on the telescope focal plane in mm.
        # detector.pixsize is in um.

        # Making wfe map...
        wfe = calc_wfe_fringe37(telescope.epd, wp['fringe37_filename'],
                                wp['reference_wl'], positions)
    else:
        print("WFE-calc method '{}' is not supported.".format(
            control_params.effect.wfe))
        exit(-1)

    return wfe


def run_calc_psf(control_params, telescope, detector, wfe):
    """Making PSF.

    Args:
        control_params: control parameters
        telescope: telescope object
        detector: detector object
        wfe: wavefront error

    Returns:
        psf: PSF in e-/sec/fp-cell for Hw=0.
    """
    # total_e_rate in e/s/m^2; wl_e_rate in um; e_rate in e/s/m^2/um.
    # these values are for an object with an apparent Hw mag of 0 mag.
    total_e_rate, wl_e_rate, e_rate = calc_response(
        control_params, telescope, detector)

    if control_params.effect.psf == 'real':
        if control_params.effect.wfe != 'fringe37':
            print('Calculating PSF...')
            psf = calc_psf(wfe, wfe.shape[0],
                           wl_e_rate, e_rate, total_e_rate,
                           telescope.total_area, telescope.aperture,
                           control_params.M_parameter, telescope.aperture.shape[0],
                           control_params.fN_parameter)
        else:
            psf = []
            for i in range(wfe.shape[0]):
                print('Calculating PSF ({}/{})...'.format(i, wfe.shape[0]))
                psf.append(calc_psf(wfe[i], wfe[i].shape[0],
                                    wl_e_rate, e_rate, total_e_rate,
                                    telescope.total_area, telescope.aperture,
                                    control_params.M_parameter, telescope.aperture.shape[0],
                                    control_params.fN_parameter))
            psf = np.array(psf)
        # psf is that of an object which has the JH color of the set value and Hw=0.
        # The unit is e/sec/fp-cell.
    elif control_params.effect.psf == 'gauss':
        print('Calculating gauss PSF...')
        psf_fwhm_arcsec = control_params.gaussPSFfwhm
        psf_fwhm_rad = psf_fwhm_arcsec*np.pi/180./3600.
        psf = calc_gauss_psf(psf_fwhm_rad, total_e_rate, telescope.total_area,
                             control_params.M_parameter, control_params.fN_parameter)
    else:
        print("The PSF mode '{}' is not supported.".format(
            control_params.effect.psf))
        exit(-1)
    return psf


def run_calc_ace(control_params, detector, ace_params):
    """Ace (Atitude control error) simulation.

    Args:
        control_params: control parameters
        detector: detector object
        ace_params: ace parameters

    Returns:
        acex: ace in x-axis, normalized with its stddev.
        acey: ace in y-axis, normalized with its stddev.
        Nts_per_plate: Number of time bins per plate.
    """
    nace = control_params.ace_control.get('nace')
    tace = control_params.ace_control.get('tace')
    print('ACE calculation mode: {}'.format(control_params.effect.ace))
    if control_params.effect.ace == 'real':
        print('  Making ACE (X)...')
        rg_acex = np.random.default_rng(
            control_params.ace_control.get('acex_seed'))
        acex, psdx = calc_ace(rg_acex, nace, tace, ace_params)
        # the standard deviation of acex is normalized to unity.

        print('  Making ACE (Y)...')
        rg_acey = np.random.default_rng(
            control_params.ace_control.get('acey_seed'))
        acey, psdy = calc_ace(rg_acey, nace, tace, ace_params)
        # the standard deviation of acey is normalized to unity.
    else:  # none/gauss mode
        print('  ACE simulation is skipped.')
        print('  Generate fake ACE(X) and ACE(Y)...')
        acex = calc_dummy_ace(np.random, nace, tace, ace_params)
        acey = calc_dummy_ace(np.random, nace, tace, ace_params)

    Nts_per_plate = int((control_params.tplate+detector.readparams.t_scan)/
                         control_params.ace_control['dtace']+0.5)
    # Number of timesteps per a plate.
    return acex, acey, Nts_per_plate


def apply_gaussian(psf, acex_std, acey_std, fp_scale):
    """In the gauss-ace mode, apply gauss filter to psf, here.

    Args:
        psf: psf in e-/sec/fp-cell.
        acex_std: std of ACE in x-axis
        acey_std: std of ACE in y-axis
        fp_scale: arcsec/fp-cell.

    Returns:
        Gaussian filtered psf
    """
    if acex_std != acey_std:
        print('In the current gauss-ace mode, acex_std must be equal to acey_std. Sorry!')
        raise ValueError('error')
    else:
        return ndimage.gaussian_filter(psf, sigma=acex_std/fp_scale)
