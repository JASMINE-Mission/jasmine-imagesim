#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""
  Make a pixel cube

  usage:
    mkpixcube.py [-h|--help] -x aX.fits -y aY.fits -v xs -w ys -p ps -n N -s tframe -f nframe -d vd [-q]

 options:
   --help       show this help message and exit
   -x aX.fits   X-axis simulated ACE file
   -v xs        aX.fits is scaled by xs
   -y aY.fits   Y-axis simulated ACE file
   -w ys        aY.fits is scaled by ys
   -n N         output image is N x N 
   -p ps        pixel scale of Output image
   -s tframe    exposure [sec] of a frame
   -f nframe    number of the frames
   -d vd        drifting velocity [pixel scale/sec]
   -q psf.fits  psf file (if not given, an analytic donuts model will be used.
""" 

from docopt import docopt             # command line interface
import numpy as np
from jis.pixsim import readflat as rf
from jis.pixsim import simpix_stable as sp
from jis.pixsim import gentraj
from jis.pixsim import makeflat as mf
from jis.jisplot import plotace 
import tqdm
import astropy.io.fits as fits
import os
import sys
import time
# Command line interface
if __name__ == '__main__':
    ts=time.time()
    args = docopt(__doc__)
    
    # Get parameters from command line
    tframe=float(args['-s'])
    nframe=int(args['-f'])
    vd=float(args['-d'])
    t=np.array(range(0,nframe))*tframe #total time sec

    x_scale = float(args['-v'])
    y_scale = float(args['-w'])
    pix_scale = float(args['-p'])
    N = int(args['-n'])

    #-----------------------------------------#
    # Loading ACE fits (should be separated in near future)
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

    if xhead["ACE-TOTT"] != yhead["ACE-TOTT"]:
        print("ACE-TOTT (total time) mismatch.")
        sys.exit(-1)
        
    Tace=xhead['ACE-TOTT']
    Nace=len(xdata)
    #-----------------------------------------#

    #artificial linear drift
    if vd>0.0:
        drift_length=vd*Tace
        #drift_azimuth=np.random.random()*2.0*np.pi
        drift_azimuth=np.pi/4.0
        drift_theta=gentraj.gentraj_drift(Nace,drift_length,drift_azimuth)
    
    #full trajectory
    theta_full=np.array([xdata*pix_scale+drift_theta[0,:],ydata*pix_scale+drift_theta[1,:]])        
    plotace.trajectory(theta_full[0,:],theta_full[1,:])
    
    Nts_per_frame= int(tframe*Nace/Tace) # number of timestep per a frame

    Nmargin=5
    Npixcube=int((np.max(theta_full)+Nmargin)*2)
    pixdim=[Npixcube,Npixcube] # adaptive pixel dimension in the aperture
    spixdim=[32,32] # subpixel dimension in a pixel

    #intrapixel
    filex="intravx.csv"
    filey="intravy.csv"
    datapath="../src/jis/pixsim/data/intrapix"
    intrapix=rf.read_intrapix(filex,filey,spixdim,datapath)

    #generate Gaussian flat 
    flat = mf.gaussian_flat(sigma=0.01)
    gpixdim=np.shape(flat) # dimension for global pixel positions

    #initial global position
    x0=(0.5*(np.shape(flat)[0]-spixdim[0]))
    y0=(0.5*(np.shape(flat)[1]-spixdim[1]))
    x=x0
    y=y0


    persistence=False
    if persistence:
        #persistence model
        tau=np.array([1.0,10.0,100.0,1000.0,10000.0]) #sec
        rho=np.array([1.e-3,1.5e-3,1.5e-3,2.e-3,3e-3])*0.8 #H2RG good detector
        NQc = len(tau) #number of the defected pixel time scale
        Qtrap = np.zeros((gpixdim[0],gpixdim[1],NQc)) #trapped charge
        xtau=tframe/tau
        Qtsave=[]
    
    lc=[]    
    jx,jy=np.int(x),np.int(y)

    pixcube=np.zeros((Npixcube,Npixcube,nframe))
    for iframe in tqdm.tqdm(range(0,nframe)):

        #the positioning system is designed to match to the future upgrade of the large drift.
        jx,jy=np.int(x),np.int(y)
        interpix=rf.flat_interpix(flat,jx,jy,pixdim,figsw=0)
        # picking temporary trajectory and local position update
        istart=iframe*Nts_per_frame
        iend=(iframe+1)*Nts_per_frame
        if iend >= Nace:
            print("insufficient time length of ACE fits.")
            sys.exit(-1)            
        theta=np.copy(theta_full[:,istart:iend])
        theta=theta+np.array([pixdim]).T/2
        pixar=sp.simpix(theta,interpix,intrapix)
        
        if persistence:
            #persistence
            Qij=read_trapped_charge(Qtrap,jx,jy,pixdim,figsw=0)        
            E0i=np.nansum(pixar,axis=2)
            Qij,Qi,Ei=persistence_const_array(xtau,rho,E0i,Qij)
            Qtrap = push_trapped_charge(Qtrap,jx,jy,pixdim,Qij)
            Qtsave.append(np.sum(Qtrap,axis=(0,1)))        
            lctmp=np.mean(np.sum(Ei))
            lc.append(lctmp)

        pixcube[:,:,iframe]=np.sum(pixar,axis=2)
            
    te=time.time()
    print(te-ts,"sec")

    #################################
    #pixcube image
    import matplotlib.pyplot as plt

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(np.sum(pixcube[:,:,:],axis=(0,1)))
    ax.set_aspect(0.7/ax.get_data_ratio())
    plt.savefig("lc.png")

    for iframe in range(0,nframe):
        plt.imshow(pixcube[:,:,iframe])
        plt.savefig("pixcube"+str(iframe)+".png")

