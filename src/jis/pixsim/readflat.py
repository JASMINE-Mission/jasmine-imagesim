import numpy as np
import pylab 
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import argparse

def flat_interpix(flat,ix,iy,pixdim,figsw=0):
    sx=ix+np.int(pixdim[0])
    sy=iy+np.int(pixdim[1])
    if(sx<np.shape(flat)[0] and sy<np.shape(flat)[1]):
        sflat=flat[ix:sx,iy:sy]
        sflat=sflat/np.mean(sflat)
    else:
        print("Invalid ix,iy")
        return None
    
    #######
    if figsw==1:
        a=plt.imshow(sflat,vmin=0.97,vmax=1.03,interpolation="nearest",cmap="bwr")
        plt.colorbar(a)
        plt.title("sigma="+str(np.std(sflat)))
        plt.savefig("flatCCD.png")
    #######

    return sflat
    


def haya2interpix(ix,iy,pixdim,dirh="/home/hirano/virtual-jasmine/pixsim/data/",file="v_flatn.fits",figsw=0):
    #ix,iy: center pixel
    hdulist = fits.open(os.path.join(dirh,file))
    flat=hdulist[0].data
    sx=ix+np.int(pixdim[0])
    sy=iy+np.int(pixdim[1])
    if(sx<np.shape(flat)[0] and sy<np.shape(flat)[1]):
        sflat=flat[ix:sx,iy:sy]
        sflat=sflat/np.mean(sflat)
    else:
        print("Invalid ix,iy")
        return None
    
    #######
    if figsw==1:
        a=plt.imshow(sflat,vmin=0.97,vmax=1.03,interpolation="nearest",cmap="bwr")
        plt.colorbar(a)
        plt.title("sigma="+str(np.std(sflat)))
        plt.savefig("flatCCD.png")
    #######

    return sflat

def read_intrapix(filex,filey,spixdim,dirh="/home/hirano/virtual-jasmine/pixsim/data/"):
    xd=np.loadtxt(os.path.join(dirh,filex),delimiter=",")
    yd=np.loadtxt(os.path.join(dirh,filey),delimiter=",")
    pxvals=np.linspace(-0.5,0.5,spixdim[0])
    xintra=np.interp(pxvals, xd.T[0],xd.T[1])
    pyvals=np.linspace(-0.5,0.5,spixdim[1])
    yintra=np.interp(pyvals, yd.T[0],yd.T[1])

    intrapix=xintra*np.array([yintra]).T
    intrapix=intrapix/np.mean(intrapix)
    return intrapix


