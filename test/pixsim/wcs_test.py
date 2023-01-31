import pytest
from jis.pixsim.wcs import pixel_scale_degree
from jis.pixsim.wcs import set_wcs
from jis.test.loadparams import load_filenames_for_test



def test_pixel_scale_radian():
    detector, control_params, telescope = load_filenames_for_test()
    ref_arcsec = 0.42441318157838753  # JASMINE pix size is approximately 0.4 arcsec
    assert pixel_scale_degree(detector,
                              telescope) * 3600 == pytest.approx(ref_arcsec)


def test_set_wcs():
    detector, control_params, telescope = load_filenames_for_test()
    jasmine_wcs = set_wcs(0.0, 0.0, 0.0, detector, telescope)
    from astropy.coordinates import SkyCoord
    coord_0 = SkyCoord('00h00m00.0s +00d00m00s', frame='fk5')
    coord = SkyCoord('00h00m00.0s +00d00m04s', frame='fk5')
    #print(coord.separation(coord_0).degree*3600) #separation is 4 arcsec
    pixels = jasmine_wcs.world_to_pixel(coord)  
    assert pixels[0] == pytest.approx(959.94604224)
    assert pixels[1] == pytest.approx(969.40333659)





if __name__ == "__main__":
    test_pixel_scale_radian()
    test_set_wcs()  