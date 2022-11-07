import numpy as np
import pytest


def load_filenames_for_test():
    import pkg_resources
    import os
    from jis.photonsim.extract_json import Detector, ControlParams, Telescope
    dirname_params = pkg_resources.resource_filename('jis', 'data/params')
    filenames = {}
    filenames['detjson'] = os.path.join(dirname_params, "det.json")
    filenames['teljson'] = os.path.join(dirname_params, "tel.json")
    filenames['ctljson'] = os.path.join(dirname_params, "ctl.json")

    detector = Detector.from_json(filenames['detjson'])
    control_params = ControlParams.from_json(filenames['ctljson'])
    telescope = Telescope.from_json(filenames['teljson'])

    return detector, control_params, telescope


def pixel_scale_radian(detector, telescope):
    """pixel scale in the unit of radian

    Args:
        detector: detector instance
        telescope: telescope instance

    Returns:
        float: pixel scale
    """
    #pixelsize is defined in the unit of micron, while efl is mm.
    return detector.pixsize * 1.e-6 / (telescope.efl * 1.e-3)


def test_pixel_scale_radian():
    detector, control_params, telescope = load_filenames_for_test()
    ref_arcsec = 0.42441318157838753  # JASMINE pix size is approximately 0.4 arcsec
    assert pixel_scale_radian(
        detector, telescope) / np.pi * 180 * 3600 == pytest.approx(ref_arcsec)


def set_wcs(ra_cen, dec_cen, detector, telescope):
    from astropy import wcs
    w = wcs.WCS(naxis=2)
    #w.wcs._naxis = [detector.npix, detector.npix] #the number of the one side pixels
    w.wcs.crpix = [detector.npix / 2.0,
                   detector.npix / 2.0]  #reference pixel coordinate
    pixel_scale = pixel_scale_radian(detector, telescope)
    w.wcs.cdelt = [pixel_scale, pixel_scale]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra_cen, dec_cen]
    
def test_set_wcs():
    import os
    filenames = load_filenames_for_test()
    detector, control_params, telescope = load_filenames_for_test()
    set_wcs(0.0,0.0, detector, telescope)
    

if __name__ == "__main__":
    test_pixel_scale_radian()
    test_set_wcs()