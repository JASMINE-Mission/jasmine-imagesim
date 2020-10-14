#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""
  Make a pixel cube

  usage:
    mkpixcube.py [-h|--help] -x aX.fits -y aY.fits -v xs -w ys -p ps -n N -o out.fits

 options:
   --help       show this help message and exit
   -x aX.fits   X-axis simulated ACE file
   -v xs        aX.fits is scaled by xs
   -y aY.fits   Y-axis simulated ACE file
   -w ys        aY.fits is scaled by ys
   -n N         output image is N x N 
   -p ps        pixel scale of Output image
   -o out.fits  output image file

""" 
from docopt import docopt             # command line interface
import numpy as np
from jis.pixsim import readflat as rf
from jis.pixsim import simpix_stable as sp
from jis.pixsim import makeflat as mf
import astropy.io.fits as fits
import os
import sys
# Command line interface
if __name__ == '__main__':
    args = docopt(__doc__)
    
    # Get parameters from command line
    x_scale = float(args['-v'])
    y_scale = float(args['-w'])
    pix_scale = float(args['-p'])
    N = int(args['-n'])
    
    xhdul = fits.open(args['-x'])
    xdata = xhdul[0].data
    xhead = xhdul[0].header
    xN = xhead['NAXIS1']
    yhdul = fits.open(args['-y'])
    ydata = yhdul[0].data
    yhead = yhdul[0].header
    yN = xhead['NAXIS1']

    if xN != yN:
        print("NAXIS1 is not the same")
        sys.exit(-1)

    pixdim=[32,32] # pixel dimension in the aperture
    spixdim=[32,32] # subpixel dimension in a pixel
    ntime=len(xdata)
    
    #intrapixel
    filex="intravx.csv"
    filey="intravy.csv"
    datapath="../src/jis/pixsim/data/intrapix"
    intrapix=rf.read_intrapix(filex,filey,spixdim,datapath)

    #generate Gaussian flat 
    flat = mf.gaussian_flat(sigma=0.01)
    gpixdim=np.shape(flat) # dimension for global pixel positions

    #trajectory
    theta=np.array([xdata,ydata])
    theta=theta*pix_scale
    nlim=10000
    theta=theta[:,0:nlim]
    print(np.shape(theta))
    nframe = 1
    tframe = 7.0 # [sec] for 1 frame
    t=np.array(range(0,nframe))*tframe #sec

    persistence=False
    if persistence:
        #persistence model
        tau=np.array([1.0,10.0,100.0,1000.0,10000.0]) #sec
        rho=np.array([1.e-3,1.5e-3,1.5e-3,2.e-3,3e-3])*0.8 #H2RG good detector
        NQc = len(tau) #number of the defected pixel time scale
        Qtrap = np.zeros((gpixdim[0],gpixdim[1],NQc)) #trapped charge
        xtau=tframe/tau
        Qtsave=[]

    #initial global position
    x0=(0.5*(np.shape(flat)[0]-spixdim[0]))
    y0=(0.5*(np.shape(flat)[1]-spixdim[1]))
    x=x0
    y=y0
    
    lc=[]    
    jx,jy=np.int(x),np.int(y)

    for iframe in range(0,nframe):
        #global position update
        #x,y=

        jxp,jyp=jx,jy
        jx,jy=np.int(x),np.int(y)
        djx,djy=jx-jxp,jy-jyp
        interpix=rf.flat_interpix(flat,jx,jy,pixdim,figsw=0)

        #local position update
        pixar=sp.simpix(theta,interpix,intrapix)

        print(np.shape(pixar))

        
        if persistence:
            #persistence
            Qij=read_trapped_charge(Qtrap,jx,jy,pixdim,figsw=0)        
            E0i=np.nansum(pixar,axis=2)
            Qij,Qi,Ei=persistence_const_array(xtau,rho,E0i,Qij)
            Qtrap = push_trapped_charge(Qtrap,jx,jy,pixdim,Qij)
            Qtsave.append(np.sum(Qtrap,axis=(0,1)))        
            lctmp=np.mean(np.sum(Ei))
            lc.append(lctmp)

