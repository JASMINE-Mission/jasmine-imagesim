import numpy as np
import argparse
import time

def addnoise(data, readnoise, rg_shot=None, rg_read=None):
    """
    Summary:
        This function adds shotnoise and readnoise.

    Args:
        data (ndarray):
            Image cube without noise in electrons.

        readnoise (float):
            Readnoise in electrons.

        rg_shot (Generator):
            Random generator for shotnoise.
            If this is not given, the shotnoise will be
            generated automatically by mk_shotnoise.

        rg_read (Generator):
            Random generator for readnoise.
            If this is not given, the readnoise will be
            generated automatically by mk_readnoise.

    Returns:
        data (ndarray)   : Image cube with noise in electrons.
        seed (dictionary): Seed created when rg_* are not given.
                           If rg_* are given, this is not returned.

    """

    seed_shot = None
    if rg_shot is None:
        shotnoise, seed_shot = mk_shotnoise(data)
    else:
        shotnoise = mk_shotnoise(data, rg_shot)

    seed_read = None
    if rg_read is None:
        readnoise, seed_read = mk_readnoise(data.shape, readnoise)
    else:
        readnoise = mk_readnoise(data.shape, readnoise, rg_read)

    data = data + shotnoise + readnoise

    if seed_shot is None and seed_read is None:
        return data
    else:
        seed = {'seed_shot': seed_shot, 'seed_read': seed_read}
        return data, seed


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

