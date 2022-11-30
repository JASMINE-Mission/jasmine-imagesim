import numpy as np


def pixel_scale_degree(detector, telescope):
    """pixel scale in the unit of degree

    Args:
        detector: detector instance
        telescope: telescope instance

    Returns:
        float: pixel scale in degree
    """
    #pixelsize is defined in the unit of micron, while efl is mm.
    return detector.pixsize * 1.e-6 / (telescope.efl * 1.e-3) / np.pi * 180


def set_wcs(ra_cen, dec_cen, rotation_angle, detector, telescope):
    """set astropy.wcs instance

    Args:
        ra_cen (float): telescope center RA
        dec_cen (float): telescope center DEC
        rotation_angle (float): rotation angle around (ra_ce, dec_cen) in degree
        detector: detector instance
        telescope: telescope instance

    Returns:
        _type_: _description_
    """
    from astropy import wcs
    w = wcs.WCS(naxis=2)
    #w.wcs._naxis = [detector.npix, detector.npix] #the number of the one side pixels
    w.wcs.crpix = [detector.npix / 2.0 + 0.5, detector.npix / 2.0 + 0.5
                   ]  #reference pixel coordinate, +1 is due to python?
    pixel_scale = pixel_scale_degree(detector, telescope)
    w.wcs.cdelt = [pixel_scale, pixel_scale]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra_cen, dec_cen]
    w.wcs.crota = [0.0,rotation_angle]
    print("pix center: ",detector.npix / 2.0, detector.npix / 2.0)
    return w
