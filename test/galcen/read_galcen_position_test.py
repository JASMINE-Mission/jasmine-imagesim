from jis.galcen.read_galcen_position import load_jscon_random_stars
from jis.galcen.read_galcen_position import random_stars_to_starplate
from jis.test.loadparams import load_filenames_for_test



def test_random_stars_to_starplate():
    import numpy as np
    from jis.pixsim.wcs import set_wcs
    detector, control_params, telescope = load_filenames_for_test()
    jasmine_wcs = set_wcs(0.0, 0.0, detector, telescope)
    random_star_data = load_jscon_random_stars()
    telescope_center_ra = np.median(random_star_data["ra"])
    telescope_center_dec = np.median(random_star_data["dec"])
    jasmine_wcs = set_wcs(telescope_center_ra, telescope_center_dec, detector,
                          telescope)
    random_star_data = load_jscon_random_stars()
    table_staplate = random_stars_to_starplate(random_star_data, jasmine_wcs)
    assert len(table_staplate) == 342157
    
def test_load_jscon_random_stars():
    data = load_jscon_random_stars()
    assert len(data["ra"]) == 342157


if __name__ == "__main__":
    test_load_jscon_random_stars()
    test_random_stars_to_starplate()