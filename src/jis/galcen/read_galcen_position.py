import pkg_resources
import pandas as pd

def load_jscon_random_stars():
    """load galcen random_stars.csv  and Hw computed by jscon (https://github.com/2ndmk2/jscon) 

    Returns:
        Pandas DataFrame
    """
    position_file = pkg_resources.resource_filename('jis',
                                                    'data/random_stars.csv')
    data = pd.read_csv(position_file, delimiter=",")
    return data

def random_stars_to_starplate(random_star_data, jasmine_wcs):
    import numpy as np
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
