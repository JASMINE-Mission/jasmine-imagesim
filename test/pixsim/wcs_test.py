import pytest
from jis.pixsim.wcs import pixel_scale_degree
from jis.pixsim.wcs import set_wcs


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


def test_pixel_scale_radian():
    detector, control_params, telescope = load_filenames_for_test()
    ref_arcsec = 0.42441318157838753  # JASMINE pix size is approximately 0.4 arcsec
    assert pixel_scale_degree(detector,
                              telescope) * 3600 == pytest.approx(ref_arcsec)


def test_set_wcs():
    filenames = load_filenames_for_test()
    detector, control_params, telescope = load_filenames_for_test()
    jasmine_wcs = set_wcs(0.0, 0.0, detector, telescope)
    from astropy.coordinates import SkyCoord
    coord = SkyCoord('00h00m00.0s +00d00m00s', frame='fk5')
    pixels = jasmine_wcs.world_to_pixel(coord)
    assert pixels[0] == pytest.approx(detector.npix / 2.)
    assert pixels[1] == pytest.approx(detector.npix / 2.)


if __name__ == "__main__":
    test_pixel_scale_radian()
    test_set_wcs()