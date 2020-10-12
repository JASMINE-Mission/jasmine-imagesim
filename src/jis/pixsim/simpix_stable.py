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

    #include "psf_ejas.h"
    #include "pixlight.h"

""",options=['-use_fast_math'])
    return source_module

def set_simpix(theta,interpix,intrapix,sigma2=2.0):

    pixdim=np.shape(interpix)
    spixdim=np.shape(intrapix)
    ntime=np.shape(theta)[1]
    thetax=theta[0,:].astype(np.float32)
    thetay=theta[1,:].astype(np.float32)
    
    Npix=pixdim[0]*pixdim[1]
    Nimg=Npix*ntime # pixel movie

    #pixlc (output pixel movie)
    pixlc=np.zeros(Nimg).astype(np.float32)
    dev_pixlc = cuda.mem_alloc(pixlc.nbytes)

    #interpix sensitivity
    finterpix=(interpix.flatten()).astype(np.float32)
    dev_interpix = cuda.mem_alloc(finterpix.nbytes)
    cuda.memcpy_htod(dev_interpix,finterpix)

    #intrapix sensitivity
    fintrapix=(intrapix.flatten()).astype(np.float32)
    dev_intrapix = cuda.mem_alloc(fintrapix.nbytes)
    cuda.memcpy_htod(dev_intrapix,fintrapix)

    dev_thetax = cuda.mem_alloc(thetax.nbytes)
    dev_thetay = cuda.mem_alloc(thetay.nbytes)
    cuda.memcpy_htod(dev_thetax,thetax)
    cuda.memcpy_htod(dev_thetay,thetay)
    
    return dev_pixlc,dev_interpix,dev_intrapix,dev_thetax,dev_thetay,pixdim,spixdim,ntime,sigma2,pixlc

def simpix(theta,interpix,intrapix,sigma2=2.0):
    start = time.time()
    #set all
    dev_pixlc,dev_interpix,dev_intrapix,dev_thetax,dev_thetay,pixdim,spixdim,ntime,sigma2,pixlc=set_simpix(theta,interpix,intrapix,sigma2)

    #kernel
    source_module = genimg()
    pkernel = source_module.get_function("pixlight")
    pkernel(dev_pixlc,dev_interpix,dev_intrapix,np.int32(ntime),\
            dev_thetax,dev_thetay,np.float32(sigma2),\
            block=(int(spixdim[0]),int(spixdim[1]),1),\
            grid=(int(pixdim[0]),int(pixdim[1])))
    cuda.memcpy_dtoh(pixlc,dev_pixlc)

    pixar=pixlc.reshape((pixdim[0],pixdim[1],ntime))
    return pixar
