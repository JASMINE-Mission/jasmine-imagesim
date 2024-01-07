import numpy as np
import json
import matplotlib.pylab as plt
from jis.pixsim import simpix_stable as sp
from jis.pixsim.addnoise import addnoise


def init_pix(filenames, control_params, detector, acex, acey, detpix_scale, driftsw):
    """Preparation for making image, Setting and plotting full trajectory.

    Args:
        filenames: names of I/O files
        control_params: control parameters
        detector: detector object
        acex: ACE in x-axis
        acey: ACE in y-axis
        detpix_scale: scale of detector pixels
        driftsw: bool if True drift included

    Returns:
        theta_full: scaled ACE in detector pixels ([[theta_x], [theta_y]]).
        pixdim: pixel dimensions of the simpix calculation.
        Npixcube: number of pixels on a side of pixcube (simpix calculation image).
    """

    if driftsw:
        from jis.photonsim.extract_json import Drift
        dft = Drift.from_json(filenames['dftjson'])
        dft.compute_drift(control_params.ace_control['dtace'], len(acex))

    # Full data of the displacement in detpix.
    # (ace[x|y] scaled and converted to detpix)
    acex_std = control_params.ace_control.get('acex_std')
    acey_std = control_params.ace_control.get('acey_std')
    if driftsw:
        theta_full = np.array([acex*acex_std/detpix_scale+dft.drift_theta[0, :],
                               acey*acey_std/detpix_scale+dft.drift_theta[1, :]])
        plt.plot(acex*acex_std/detpix_scale +
                 dft.drift_theta[0, :], acey*acey_std/detpix_scale+dft.drift_theta[1, :], '.')
        plt.savefig(filenames['thetapng'])
    else:
        theta_full = np.array(
            [acex*acex_std/detpix_scale, acey*acey_std/detpix_scale])

    Npixcube = int((np.max(np.abs(theta_full))+detector.nmargin)*2)
    pixdim = [Npixcube, Npixcube]  # adaptive pixel dimension in the aperture.
    return theta_full, pixdim, Npixcube


def uniform_flat(detector):
    """make uniform inter/intra pixel.

    Args:
        detector: detector object

    Returns:
        uniform_flat_interpix
        uniform_flat_intrapix
    """
    uniform_flat_interpix = np.ones_like(detector.flat.interpix)
    uniform_flat_intrapix = np.ones_like(detector.flat.intrapix)
    return uniform_flat_interpix, uniform_flat_intrapix


def init_images(control_params, detector, prior_dark = True, addnoise=True):
    """initialize pixcube.

    Args:
        control_params: control parameters
        detector: detector object
        prior_dark: if the dark is added (True) or not (False). default: True
        addnoise: switch for noise-addition function.

    Returns:
        global pixel cube images
    """
    if prior_dark:
        pixcube_global = global_dark(control_params, detector, addnoise=addnoise)
    else:
        pixcube_global = np.zeros(shape=(detector.npix, detector.npix, control_params.nplate))
    
    return pixcube_global

def global_dark(control_params, detector, addnoise=True, digitize=True):
    """compute global dark image

    Args: 
        control_params: control parameters
        detector: detector object
        addnoise: switch for noise-addition function.
        digitize: switch for digitization function.

    Returns:
        global pixel cube dark image.
        if digitize=True, the unit is adu.
        if digitize=False, the unit is e-.
    """
    pixcube_global_dark = np.zeros(shape=(detector.npix, detector.npix, control_params.nplate))
    pixcube_global_dark += detector.idark * control_params.tplate
    if addnoise:
        pixcube_global_dark, seed = addnoise(pixcube_global_dark, np.sqrt(2.)*detector.readnoise)
    # Digitization: converting to adu/pix/plate.
    if digitize:
        pixcube_global_dark = np.round(pixcube_global_dark/detector.gain)
    return pixcube_global_dark


def set_positions(line, Npixcube):
    """Position setting.

    Args:
        line: each of table_starplate
        Npixcube: number of one side of pixcube

    Returns:
        xc_local, yc_local, x0_global, y0_global, xc_global, yc_global
    """
    xc_global = line['x pixel'] - 1  # Stellar pos. in glob. coord (X).
    yc_global = line['y pixel'] - 1  # Stellar pos. in glob. coord (Y).
    # Origin pix position in global coord (x).
    x0_global = int(xc_global - Npixcube*0.5 + 0.5)
    # Origin pix position in global coord (y).
    y0_global = int(yc_global - Npixcube*0.5 + 0.5)
    xc_local = xc_global - x0_global  # Stellar position (local; x).
    yc_local = yc_global - y0_global  # Stellar position (local; y).
    return xc_local, yc_local, x0_global, y0_global, xc_global, yc_global


def make_local_flat(control_params, detector, x0_global, y0_global, pixdim):
    """Making local flat data.

    Args:
        control_params: control parameters
        detector: detector object
        x0_global: x global position
        y0_global: y global position
        pixdim: pixel dimension

    Returns:
        local interpixel fluctuation
    """
    from jis.pixsim import readflat as rf
    if control_params.effect.flat_interpix is True:
        interpix_local = rf.flat_interpix(
            detector.flat.interpix, x0_global, y0_global, pixdim, figsw=0)
    else:
        uniform_flat_interpix, uniform_flat_intrapix = uniform_flat(detector)
        interpix_local = rf.flat_interpix(
            uniform_flat_interpix, x0_global, y0_global, pixdim, figsw=0)
    return interpix_local


def index_control_trajectory(control_params, iplate, Nts_per_plate):
    """get indices for trajectory.

    Args:
        control_params: control parameters
        iplate: index for plate
        Nts_per_plate: number of time bins per plante

    Returns:
        istart for trajectory
        iend for trajectory
    """
    istart = iplate * Nts_per_plate
    iend = (iplate+1)*Nts_per_plate
    # In no-ace mode, we make a single image with simpix to reduce the calculation time.
    # Below is a trick for that. After executing simpix, we will copy it to make Nts_per_plate shots.
    if control_params.effect.ace != 'real':
        iend = istart+1
    return istart, iend


def calc_theta(theta_full, istart, iend, xc_local, yc_local):
    """calc theta (star central position)

    Args:
        theta_full: full trajectory
        istart: start index by index_control_trajectory
        iend: end index by index_control_trajectory
        xc_local: x center position for local
        yc_local: y center position for local

    Returns:
       theta
    """

    # Displacement from the initial position.
    theta = np.copy(theta_full[:, istart:iend])
    # Displacement in local coord.
    theta = theta + np.array([[xc_local, yc_local]]).T
    # 0.5-pix shift to treat the coodinate difference.
    theta = theta + np.array([[0.5, 0.5]]).T
    # Global coord.   : (0, 0) is the center of the bottom-left corner pixel.
    # Local coord.    : (0, 0) is the center of the bottom-left corner pixel..
    # Coord. in simpix: (0, 0) is the bottom-left corner of the bottom-left corner pixel.
    return theta


def run_simpix(control_params, theta, interpix_local, flat_intrapix, psfarr, psfcenter, psfscale, Nts_per_plate):
    """run simpix and normalize it

    Args:
        control_params: control parameters       
        theta: trajectory
        interpix_local: (local) interpixel fluctuation
        flat_intrapix: intrapix fluctuation
        psfarr: psfarr
        psfcenter: psf center
        psfscale: psfscale
        Nts_per_plate: number of time bins per plate


    Returns:
        unscaled pixar
    """

    pixar = sp.simpix(theta, interpix_local, flat_intrapix,
                      psfarr=psfarr, psfcenter=psfcenter, psfscale=psfscale)\
                      / (psfscale*psfscale)*control_params.ace_control['dtace']/(1./Nts_per_plate)

    # pixar is in e/pix/control_params.ace_control['dtace'].
    # In none/gauss mode, we copy the single-shot image to make the full-movie cube.
    if control_params.effect.ace != 'real':
        upixar = pixar[:, :, 0]
        nxt, nyt = np.shape(upixar)
        pixar = upixar[:, :, np.newaxis]+np.zeros((nxt, nyt, Nts_per_plate))
        pixar = pixar/Nts_per_plate
    return pixar
        
def scaling_pixar(pixar, mag):
    """scale pixar.

    Args:
        pixar: pixar
        mag: magnitude

    Returns:
        pixar
    """

    pixar = pixar * 10.**(mag/(-2.5))  # magnitude scaling.

    return pixar


def add_varability(pixar, injlc_iplate):
    """add variability to pixar.

    Args:
        pixar: pixel array
        injlc_iplate: injlc[iplate]

    Returns:
        pixel array variability imprinted
    """
    return pixar*injlc_iplate


def add_dark_current(control_params, detector, pixar):
    """Adding dark current (including stray light).

    Args:
        control_params: control parameters
        detector: detector object
        pixar:  pixel array

    Returns:
        pixel array dark current added
    """
    dark = np.ones(shape=pixar.shape) * detector.idark * \
        control_params.ace_control['dtace']
    return pixar + dark
