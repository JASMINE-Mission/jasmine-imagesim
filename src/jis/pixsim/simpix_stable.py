import pycuda.autoinit
import pycuda.driver as cuda
import pycuda.compiler
from pycuda.compiler import SourceModule
import time
import numpy as np

def genimg():
    source_module = SourceModule("""
    #define NMXCACHE 1024 
    #define PI 3.14159265359

    __shared__ float cache[NMXCACHE];    

    #include "psf_donut.h"
    #include "pixlight.h"

""",options=['-use_fast_math'])
    return source_module

def set_simpix(theta,interpix,intrapix,sigma2=2.0):
    """
    Summary:
        This function makes preparations for simpix.
        (Memory allocation etc.)

    Args:
        theta    (ndarray): Time-series data of the PSF center position
                            (pix; [0,:] is x_center; [1,:] is y_center).
        interpix (ndarray): Interpixel map (2-d array).
        intrapix (ndarray): Intrapixel map (2-d array).
        sigma2   (float)  : sigma^2 of gaussian PSF (pix^2).
                            Default: 2.0

    Return:
        dev_pixlc    (DeviceAllocation): Memory for output movie data (pixlc).
        dev_interpix (DeviceAllocation): Memory for interpix pattern map.
        dev_intrapix (DeviceAllocation): Memory for intrapix pattern map.
        dev_thetax   (DeviceAllocation): Memory for psf-center in x-dir. (thetax).
        dev_thetay   (DeviceAllocation): Memory for psf-center in y-dir. (thetay).
        pixdim       (ndarray)         : Dimensions of interpix.
        spixdim      (ndarray)         : Dimensions of intrapix.
        ntime        (int)             : Number of time grid points.
        sigma2       (float)           : sigma^2 of gaussian PSF (pix^2).
        pixlc        (ndarray)         : Array for output movie data.
 
    """

    pixdim  = np.shape(interpix)
    spixdim = np.shape(intrapix)
    ntime   = np.shape(theta)[1]
    thetax  = theta[0,:].astype(np.float32)
    thetay  = theta[1,:].astype(np.float32)
    
    Npix = pixdim[0] * pixdim[1]
    Nimg = Npix * ntime # pixel movie

    # pixlc (output pixel movie)
    pixlc = np.zeros(Nimg).astype(np.float32)
    dev_pixlc = cuda.mem_alloc(pixlc.nbytes)

    # interpix sensitivity
    finterpix = (interpix.flatten()).astype(np.float32)
    dev_interpix = cuda.mem_alloc(finterpix.nbytes)
    cuda.memcpy_htod(dev_interpix,finterpix)

    # intrapix sensitivity
    fintrapix = (intrapix.flatten()).astype(np.float32)
    dev_intrapix = cuda.mem_alloc(fintrapix.nbytes)
    cuda.memcpy_htod(dev_intrapix,fintrapix)

    dev_thetax = cuda.mem_alloc(thetax.nbytes)
    dev_thetay = cuda.mem_alloc(thetay.nbytes)
    cuda.memcpy_htod(dev_thetax,thetax)
    cuda.memcpy_htod(dev_thetay,thetay)
    
    return dev_pixlc, dev_interpix, dev_intrapix, dev_thetax, dev_thetay,\
           pixdim, spixdim, ntime, sigma2, pixlc


def simpix(theta, interpix, intrapix, sigma2=2.0, readnoise=15.):
    """
    Summary:
        This function makes a movie data 
        based on the time-series data of the observed position (theta)
        with taking interpix/intrapix flat data into account.
        Calculation is done by pixlight command.
        The PSF is assumed to be a gaussian function.

    Args:
        theta     (ndarray): Time-series data of the PSF position.
        interpix  (ndarray): Interpixel flat pattern (2-d array).
        intrapix  (ndarray): Intrapixel flat pattern (2-d array).
        sigma2    (float)  : sigma^2 of gaussian PSF.
        readnoise (float)  : Readnoise in electrons (default 15e-).

    Returns:
        pixar (ndarray): Calculated movie data (3-d array).  

    """

    # sigma2 is dummy
    start = time.time()
    # set all
    dev_pixlc, dev_interpix, dev_intrapix, dev_thetax, dev_thetay,\
    pixdim, spixdim, ntime, sigma2, pixlc = set_simpix(theta, interpix, intrapix, sigma2)

    #kernel
    source_module = genimg()
    pkernel = source_module.get_function("pixlight")
    pkernel(dev_pixlc, dev_interpix, dev_intrapix, np.int32(ntime),\
            dev_thetax, dev_thetay, np.float32(sigma2),\
            block=(int(spixdim[0]), int(spixdim[1]),1),\
            grid=(int(pixdim[0]),int(pixdim[1])))

    cuda.memcpy_dtoh(pixlc,dev_pixlc)

    pixar = pixlc.reshape((pixdim[0], pixdim[1], ntime))

    return pixar
