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

