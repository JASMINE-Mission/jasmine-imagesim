from jis.galcen.read_galcen_position import load_jscon_random_stars


def radec_to_starpalte(random_star_data, jasmine_wcs):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    import pandas as pd
    coord = SkyCoord(ra=random_star_data["ra"].values * u.degree,
                     dec=random_star_data["dec"].values * u.degree)
    xpixel, ypixel = jasmine_wcs.world_to_pixel(coord)
    hwmag = random_star_data["hwmag"].values
    lambda_ = np.zeros_like(hwmag)
    beta_ = np.zeros_like(hwmag)
    plate_index = np.zeros_like(hwmag)
    star_index = np.array(range(len(hwmag)))
    table_starplate = pd.DataFrame({
        "star index": star_index,
        "x pixel": xpixel,
        "y pixel": ypixel,
        "lambda": lambda_,
        "beta": beta_,
        "Hwmag": hwmag
    })
    return table_starplate

def test_radec_to_starpalte():
    import numpy as np
    from jis.pixsim.wcs import set_wcs
    telescope_center_ra = np.median(random_star_data["ra"])
    telescope_center_dec = np.median(random_star_data["dec"])
    jasmine_wcs = set_wcs(telescope_center_ra, telescope_center_dec, detector,
                          telescope)
    random_star_data = load_jscon_random_stars()
    radec_to_starpalte(random_star_data, jasmine_wcs)
    
    
def test_load_jscon_random_stars():
    data = load_jscon_random_stars()
    assert len(data["ra"]) == 342157


if __name__ == "__main__":
    test_load_jscon_random_stars()