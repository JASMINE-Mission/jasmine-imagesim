import pkg_resources
import numpy as np


def load_jscon_data():
    """load galcen star positions and Hw computed by jscon (https://github.com/2ndmk2/jscon) 

    Returns:
        _type_: _description_
    """
    position_file = pkg_resources.resource_filename('jis',
                                                    'data/jscon_position.npz')
    gal_l, gal_b, hw, hwtarget = np.load(position_file)["arr_0"]
    return gal_l, gal_b, hw, hwtarget
