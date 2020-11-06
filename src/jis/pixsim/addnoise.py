import numpy as np
import argparse
import time

def mk_shotnoise(data, rg=None):
    """
    Summary:
        This function returns shotnoise cube.

    Args:
        data (ndarray):
             Image cube in electrons.

        rg (Generator):
             Random generator. If this is not given, 
             the noise will be generated by using np.random.PCG64
             with a seed created automatically from time.

    Returns:
        shotnoise (ndarray):
             Calculated shotnoise cube in electrons.

        seed (int):
             Seed created when rg is not given.
             If rg is given, this is not returned.

    """

    seed = None
    if rg is None:
        seed = round(time.time()*10)
        rg   = np.random.Generator(np.random.PCG64(seed))

    sigma     = np.sqrt(data) 
    shotnoise = rg.standard_normal(data.shape)*sigma

    if seed is None:
        return shotnoise
    else:
        return shotnoise, seed


def mk_readnoise(shape, sigma, rg=None):
    """
    Summary:
        This function returns readnoise cube.

    Args:
        shape (array_like):
            Shape of the cube.

        sigma (float):
            One-sigma of the readnoise in electrons.

        rg (Generator):
            Random generator. If this is not given,
            the noise will be generated by using np.random.PCG64
            with a seed created automatically from time.

    Returns:
        readnoise (ndarray):
            Calculated readnoise cube in electrons.

        seed (int):
            Seed created when rg is not given.
            If rg is given, this is not returned.
    """

    seed = None
    if rg is None:
        seed = round(time.time()*10)
        rg   = np.random.Generator(np.random.PCG64(seed))

    readnoise = rg.standard_normal(shape)*sigma

    if seed is None:
        return readnoise
    else:
        return readnoise, seed
    

### Obsolete? (TK) ######################################    
#def gausian_ppm(arshape,ppm=500.0):
#    sigma=ppm*1.e-6
#    print(arshape)
#    noise=np.random.standard_normal(arshape)*sigma
#    return noise

def ppm2nphoton(ppm):
    return 1.0/(ppm*1.e-6)**2

def convfactor_photon(ppm,meanar):
    cfac=ppm2nphoton(ppm)/meanar
    return cfac

