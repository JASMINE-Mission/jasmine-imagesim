import pycuda.autoinit
import pycuda.driver as cuda
import pycuda.compiler
from pycuda.compiler import SourceModule
import time
import numpy as np

def genimg_donut(spixdim):
    print("Analytic donut PSF model used.")
    source_module = SourceModule(
    "#define NMXCACHE "+str(spixdim[0]*spixdim[1])+"\n"+ 
    """
    #define PI 3.14159265359
    
    __shared__ float cache[NMXCACHE];    
    
    #include "psf_donut.h"
    #include "pixlight_analytic.h"
    
    """,options=['-use_fast_math'])
    
    return source_module

def genimg_custom(spixdim,Nsubtilex, Nsubtiley):        
    print("Custom PSF model used.")
    source_module = SourceModule(
    "#define NMXCACHE "+str(spixdim[0]*spixdim[1]+Nsubtilex*Nsubtiley)+"\n"+ 
    """
    #define NMXCACHE 1024 
    #define PI 3.14159265359
    
    __shared__ float cache[NMXCACHE];    
    
    #include "psf_custom.h"
    #include "pixlight_custom.h"
    
    """,options=['-use_fast_math'])
        
    return source_module

def set_simpix(theta,interpix,intrapix,sigma2=2.0):
    """
    Summary:
        This function makes preparations for simpix.
        (Gloabl memory allocation etc.)

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

def pix2psfpix(pixpos,theta,psfcenter,psfscale):
    # psfpos [psf_pix]: psfarr pixel position as a function of (detector) pixel position (pixpos [pix])
    # theta: the PSF center position [pix]
    # psfscale [pix/psf_pix]
    
    psfpos=psfcenter + (pixpos - theta)/psfscale
    return psfpos

def set_custom(theta,psfarr,psfcenter,psfside,pixdim):
    # psf array
    fpsfarr = (psfarr.flatten()).astype(np.float32)
    dev_psfarr = cuda.mem_alloc(fpsfarr.nbytes)
    cuda.memcpy_htod(dev_psfarr,psfarr)

    #pix and psfarr center
    pixcenter=pixdim/2.0
    psfdim=np.shape(psfarr)
    psfscale=psfside/psfdim #[pix/psf_pix] 1 psfarr pixel scale in the unit of (detector) pixel scale 
    
    subtilex=np.zeros(pixdim)
    subtiley=np.zeros(pixdim)

    thetamax=np.max(theta,axis=1)
    thetamin=np.min(theta,axis=1)
    
    Nsubtilex=0
    Nsubtiley=0

    for ix in range(0,pixdim[0]):
        for iy in range(0,pixdim[1]):
            minpixpos=thetamin+np.array([ix,iy])
            minpsfpos=pix2psfpix(minpixpos,theta,psfcenter,psfscale)
            subtilex[ix,iy]=minpsfpos[0]-1.0/psfscale[0]/2-1
            subtiley[ix,iy]=minpsfpos[1]-1.0/psfscale[1]/2-1

            maxpixpos=thetamax+np.array([ix,iy])
            maxpsfpos=pix2psfpix(maxpixpos,theta,psfcenter,psfscale)

            #check the maximum size of subtile
            smax=maxpsfpos[0]+1.0/psfscale/2+1        
            Nsx=smax[0]-subtilex[ix,iy]
            if Nsx > Nsubtilex:
                Nsubtilex=Nsx                
            Nsy=smax[1]-subtiley[ix,iy]
            if Nsy > Nsubtiley:
                Nsubtiley=Nsy

                
    subtilex = (np.array(subtilex)).astype(np.float32)
    dev_subtilex = cuda.mem_alloc(subtilex.nbytes)
    cuda.memcpy_htod(dev_subtilex,subtilex)

    subtiley = (np.array(subtiley)).astype(np.float32)
    dev_subtiley = cuda.mem_alloc(subtiley.nbytes)
    cuda.memcpy_htod(dev_subtiley,subtiley)
            
    return dev_psfarr, dev_subtilex, dev_subtiley, Nsubtilex, Nsubtiley


def simpix(theta, interpix, intrapix, sigma2=2.0, psfarr=None, psfcenter=None, psfside=None):
    """
    Summary:
        This function makes a movie data 
        based on the time-series data of the observed position (theta)
        with taking interpix/intrapix flat data into account.
        Calculation is done by pixlight command.
        The PSF is assumed to be a gaussian function.

    Args:
        theta    (ndarray): Time-series data of the PSF position.
        interpix (ndarray): Interpixel flat pattern (2-d array).
        intrapix (ndarray): Intrapixel flat pattern (2-d array).
        sigma2   (float)  : sigma^2 of gaussian PSF.
        psf      (ndarray): psf array or None=the analytic donut.
        psfside  (float)  : psf array size in the unit of (detector) pixel.
    Returns:
        pixar (ndarray): Calculated movie data (3-d array).  

    """

    start = time.time()
    
    # set all
    dev_pixlc, dev_interpix, dev_intrapix, dev_thetax, dev_thetay,\
        pixdim, spixdim, ntime, sigma2, pixlc = set_simpix(theta, interpix, intrapix, sigma2)     # sigma2 is dummy

    #kernel
    if psf is None:
        source_module = genimg_donut(spixdim)
        pkernel = source_module.get_function("pixlight_analytic")
        pkernel(dev_pixlc, dev_interpix, dev_intrapix, np.int32(ntime),\
                dev_thetax, dev_thetay, np.float32(sigma2),\
                block=(int(spixdim[0]), int(spixdim[1]),1),\
                grid=(int(pixdim[0]),int(pixdim[1])))
    else:
        dev_psfarr, dev_subtilex, dev_subtiley, Nsubtilex, Nsubtiley = set_custom(theta, psfarr, psfcenter, psfside, pixdim)
        source_module = genimg_custom(spixdim,Nsubtilex, Nsubtiley)
        pkernel = source_module.get_function("pixlight_custom")
        pkernel(dev_pixlc, dev_interpix, dev_intrapix, dev_psfarr,\
                np.int32(ntime),np.int32(Nsubtilex),np.int32(Nsubtiley),\
                dev_thetax, dev_thetay,\
                block=(int(spixdim[0]), int(spixdim[1]),1),\
                grid=(int(pixdim[0]),int(pixdim[1])))

        
    cuda.memcpy_dtoh(pixlc,dev_pixlc)

    pixar = pixlc.reshape((pixdim[0], pixdim[1], ntime))

    return pixar

