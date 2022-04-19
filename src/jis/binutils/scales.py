import numpy as np


def get_pixelscales(control_params, telescope, detector):
    """ 

    Args:
        control_params: control parameters
        telescope: telescope object
        detector: detector object

    Returns:
        det. pix. scale in arcsec/pix
        in rad/fp-cell.   
        arcsec/fp-cell. 
        det-pix/fp-cell.

    """
    detpix_scale = detector.pixsize*1.e-6/telescope.efl/1.e-3*180.*3600./np.pi
    fp_cellsize_rad = (1./control_params.M_parameter)*1.e-3
    fp_scale = fp_cellsize_rad * 3600.*180./np.pi
    psfscale = fp_scale/detpix_scale
    return detpix_scale, fp_cellsize_rad, fp_scale, psfscale


def get_tday(control_params, detector):
    """t day.

    Args:
       control_params: control parameters
        detector: detector object

    Returns:
        tday
    """
    return (control_params.tplate+detector.readparams.t_scan)*np.array(range(0, control_params.nplate))/3600/24
