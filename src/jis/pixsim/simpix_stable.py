import pycuda.autoinit
import pycuda.driver as cuda
import pycuda.compiler
from pycuda.compiler import SourceModule
import time
import numpy as np
def genimg_donut(spixdim):
    cudacode=\
    "    #define NMXCACHE "+str(spixdim[0]*spixdim[1])+"\n"\
    +"""
    #define PI 3.14159265359
    
    __shared__ float cache[NMXCACHE];    
    
    #include "psf_donut.h"
    #include "pixlight_analytic.h"
    
    """
    source_module = SourceModule(cudacode,options=['-use_fast_math'])
    
    return source_module

def genimg_custom(psfdim, spixdim, psfcenter, psfscale, Nsubtilex, Nsubtiley):        
    cudacode=\
    "    #define NMXCACHE "+str(spixdim[0]*spixdim[1]+Nsubtilex*Nsubtiley)+"\n"\
    +"    #define NNSUBTILE "+str(Nsubtilex*Nsubtiley)+"\n"\
    +"    #define NINTRA "+str(spixdim[0]*spixdim[1])+"\n"\
    +"    #define NSUBTILEX "+str(Nsubtiley)+"\n"\
    +"    #define NSUBTILEY "+str(Nsubtilex)+"\n"\
    +"    #define PSFDIMX "+str(psfdim[1])+"\n"\
    +"    #define PSFDIMY "+str(psfdim[0])+"\n"\
    +"    #define PSFCENTERX "+str(psfcenter[1])+"\n"\
    +"    #define PSFCENTERY "+str(psfcenter[0])+"\n"\
    +"    #define PSFSCALE "+str(psfscale)+"\n"\
    +"""
    #define PI 3.14159265359
    
    __shared__ float cache[NMXCACHE];    
    
    #include "pixlight_custom.h"
    
    """
    source_module = SourceModule(cudacode,options=['-use_fast_math'])

    return source_module

def set_simpix(theta,interpix,intrapix):
    """
    Summary:
        This function makes preparations for simpix.
        (Gloabl memory allocation etc.)

    Args:
        theta    (ndarray): Time-series data of the PSF center position
                            (pix; [0,:] is x_center; [1,:] is y_center).
        interpix (ndarray): Interpixel map (2-d array).
        intrapix (ndarray): Intrapixel map (2-d array).

    Return:
        dev_pixlc    (DeviceAllocation): Memory for output movie data (pixlc).
        dev_interpix (DeviceAllocation): Memory for interpix pattern map.
        dev_intrapix (DeviceAllocation): Memory for intrapix pattern map.
        dev_thetax   (DeviceAllocation): Memory for psf-center in x-dir. (thetax).
        dev_thetay   (DeviceAllocation): Memory for psf-center in y-dir. (thetay).
        pixdim       (ndarray)         : Dimensions of interpix.
        spixdim      (ndarray)         : Dimensions of intrapix.
        ntime        (int)             : Number of time grid points.
        pixlc        (ndarray)         : Array for output movie data.
 
    """

    pixdim  = np.array(np.shape(interpix))
    spixdim = np.array(np.shape(intrapix))
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
           pixdim, spixdim, ntime, pixlc

def pix2psfpix(pixpos,pixcenter,psfcenter,psfscale):
    # psfpos [fp-cell]: psfarr pixel position as a function of (detector) pixel position (pixpos [pix])
    # pixcenter: the detector center position [pix,pix]
    # psfscale [pix/fp-cell]
    
    psfpos=psfcenter + (pixpos - pixcenter)/psfscale
    return psfpos

def set_custom(theta,psfarr,psfcenter,psfscale,pixdim,spixdim):
    import sys

    # psf array
    fpsfarr = (psfarr.flatten()).astype(np.float32)
    dev_psfarr = cuda.mem_alloc(fpsfarr.nbytes)
    cuda.memcpy_htod(dev_psfarr,fpsfarr)

    #subtile initialization
    subtilex=np.zeros(pixdim)
    subtiley=np.zeros(pixdim)

    thetamax=np.max(theta,axis=1)
    thetamin=np.min(theta,axis=1)
    thetamed=np.median(theta,axis=1)

    # psfscale [pix/fp-cell]
    full_pixsize=1.0/psfscale
    psfdim=np.array(np.shape(psfarr))
    Nsubtilex=0
    Nsubtiley=0

    #checking shared memory size    
    dev=cuda.Device(0)
    shared_memory_size=dev.MAX_SHARED_MEMORY_PER_BLOCK
    
    for ix in range(0,pixdim[0]):
        for iy in range(0,pixdim[1]):
            pixpos=np.array([ix,iy])
            maxpsfpos=pix2psfpix(pixpos,thetamin,psfcenter,psfscale)
            minpsfpos=pix2psfpix(pixpos,thetamax,psfcenter,psfscale)

            jx=int(minpsfpos[0]-full_pixsize-1)
            if jx >= 0 and jx < psfdim[0]:
                subtilex[ix,iy]=jx
            else:
                sys.exit("OVER jx")

            jy=int(minpsfpos[1]-full_pixsize-1)
            if jy >= 0 and jy < psfdim[1]:
                subtiley[ix,iy]=jy
            else:
                sys.exit("OVER jy")
            
            #check the maximum size of subtile
            smax=maxpsfpos+full_pixsize+1
            Nsx=smax[0]-subtilex[ix,iy]+1
            if Nsx > Nsubtilex:
                Nsubtilex=int(Nsx)                
            Nsy=smax[1]-subtiley[ix,iy]+1
            if Nsy > Nsubtiley:
                Nsubtiley=int(Nsy)
                
    subtilex_flat = (np.array(subtilex)).astype(np.int32)    
    dev_subtilex = cuda.mem_alloc(subtilex_flat.nbytes)
    cuda.memcpy_htod(dev_subtilex,subtilex_flat)

    subtiley_flat = (np.array(subtiley)).astype(np.int32)
    dev_subtiley = cuda.mem_alloc(subtiley_flat.nbytes)
    cuda.memcpy_htod(dev_subtiley,subtiley_flat)

    if (Nsx*Nsy+spixdim[0]*spixdim[1]) > shared_memory_size:
        print("your system has "+str(shared_memory_size)+" byte/block for shared memory.")
        print("But we need "+str(Nsx*Nsy+spixdim[0]*spixdim[1])+" byte/block.")
        sys.exit("Error: Shared memory size not enough.")

    
    return dev_psfarr, dev_subtilex, dev_subtiley, subtilex, subtiley, Nsubtilex, Nsubtiley, psfdim


def emurate_pixlight_custom(theta_instant,psfarr,pixdim,spixdim,subtilex,subtiley,Nsubtilex, Nsubtiley,psfcenter,psfscale):
    #------------------------------------------------------
    ## emurating cuda pixlight_custom.h
    testpsf=np.zeros(pixdim)
    allpsf=np.zeros(pixdim*spixdim)
    
    for ix in range(0,pixdim[0]):
        for iy in range(0,pixdim[1]):
            jx=int(subtilex[ix,iy])
            jy=int(subtiley[ix,iy])
#            subtile=psfarr[jx:jx+Nsubtilex,jy:jy+Nsubtiley]
            subtile=psfarr[jx:jx+Nsubtilex,jy:jy+Nsubtiley]
            if jx<0 or jy<0:
                testpsf[ix,iy]=None
                print("out of range",ix,iy)
            else:
                testpsf[ix,iy]=psfarr[jx,jy]
                for kx in range(0,spixdim[0]):
                    for ky in range(0,spixdim[1]):

                        pixpos=np.array([ix+kx/spixdim[0],iy+ky/spixdim[1]])
                        psfpos=pix2psfpix(pixpos,theta_instant,psfcenter,psfscale)
                        
                        #bilinear interpolation
                        x=psfpos[0]-jx #x-position in subtile
                        y=psfpos[1]-jy #y-position in subtile
                        x1=int(x)
                        x2=x1+1
                        y1=int(y)
                        y2=y1+1

#                        Q11=psfarr[x1+jx,y1+jy]
#                        Q12=psfarr[x1+jx,y2+jy]
#                        Q21=psfarr[x2+jx,y1+jy]
#                        Q22=psfarr[x2+jx,y2+jy]
                        
                        Q11=subtile[x1,y1]
                        Q12=subtile[x1,y2]
                        Q21=subtile[x2,y1]
                        Q22=subtile[x2,y2]
                        
                        F1=(x2-x)*Q11 + (x-x1)*Q21
                        F2=(x2-x)*Q12 + (x-x1)*Q22

                        lx=spixdim[0]*ix + kx
                        ly=spixdim[1]*iy + ky
                        allpsf[lx,ly]=(y2-y)*F1 + (y-y1)*F2
                            
    import matplotlib.pyplot as plt
    #a=plt.imshow(testpsf)
    a=plt.imshow(allpsf,interpolation=None,cmap="rainbow")
    plt.colorbar(a)
    plt.show()
    sys.exit()
    return allpsf
    #------------------------------------------------------




def simpix(theta, interpix, intrapix, psfarr=None, psfcenter=None, psfscale=None):
  
    """
    Summary:
        This function makes a movie data 
        based on the time-series data of the observed position (theta)
        with taking interpix/intrapix flat data into account.
        Calculation is done by pixlight command.
        If psfarr,psfcenter,psfscale are given, the bilinear interporation of psfarr provides the psf. If not, the analytic psf assuming a donut shape will be applied.

    Args:
        theta     (ndarray): Time-series data of the PSF position.
        interpix  (ndarray): Interpixel flat pattern (2-d array).
        intrapix  (ndarray): Intrapixel flat pattern (2-d array).
        psfarr    (ndarray): psf array or None=the analytic donut.
        psfcenter (ndarray): psf center in the unit of fp-cell
        psfscale  (float)  : psf pixel-size in the unit of (detector) pix [pix/fp-cell]
        readnoise (float)  : Readnoise in electrons (default 15e-).

    Returns:
        pixar (ndarray): Calculated movie data (3-d array).  

    """

    start = time.time()
    
    # set all
    dev_pixlc, dev_interpix, dev_intrapix, dev_thetax, dev_thetay,\
        pixdim, spixdim, ntime, pixlc = set_simpix(theta, interpix, intrapix) #
    #kernel
    if psfarr is None:
        source_module = genimg_donut(spixdim)
        pkernel = source_module.get_function("pixlight_analytic")
        pkernel(dev_pixlc, dev_interpix, dev_intrapix, np.int32(ntime),\
                dev_thetax, dev_thetay,\
                block=(int(spixdim[0]), int(spixdim[1]),1),\
                grid=(int(pixdim[0]),int(pixdim[1])))
    else:
        if psfscale is None or psfcenter is None:
            import sys
            sys.exit("Error: Provide both psfscale and psfcenter in simpix.")

        dev_psfarr,dev_subtilex, dev_subtiley, subtilex, subtiley, Nsubtilex, Nsubtiley, psfdim =\
        set_custom(theta, psfarr, psfcenter, psfscale, pixdim, spixdim)
#        emurate_pixlight_custom(theta[:,0],psfarr,pixdim,spixdim,subtilex,subtiley,Nsubtilex, Nsubtiley,psfcenter,psfscale)
        
        source_module = genimg_custom(psfdim, spixdim, psfcenter, psfscale, Nsubtilex, Nsubtiley)        
        pkernel = source_module.get_function("pixlight_custom")
        pkernel(dev_pixlc, dev_interpix, dev_intrapix,\
                dev_psfarr,dev_subtilex, dev_subtiley,\
                np.int32(ntime),dev_thetax, dev_thetay,\
                block=(int(spixdim[0]), int(spixdim[1]),1),\
                grid=(int(pixdim[0]),int(pixdim[1])))

        
    cuda.memcpy_dtoh(pixlc,dev_pixlc)

    pixar = pixlc.reshape((pixdim[0], pixdim[1], ntime))

    return pixar

