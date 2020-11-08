#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""
  Make a pixel cube

  usage:
    mkpixcube.py [-h|--help] [-l lc.fits] -x aX.fits -y aY.fits -v xs -w ys -p ps -n N -s tframe -f nframe -d vd --det det.json -o pixcube [-m] [--persistence] [--psf] psf.fits 

 options:
   --help         show this help message and exit
   -l lc.fits     light curve
   -x aX.fits     X-axis simulated ACE file
   -v xs          aX.fits is scaled by xs
   -y aY.fits     Y-axis simulated ACE file
   -w ys          aY.fits is scaled by ys
   -n N           output image is N x N
   -p ps          pixel scale of Output image
   -s tframe      exposure [sec] of a frame
   -f nframe      number of the frames
   -d vd          drifting velocity [pixel scale/sec]
   --det det.json detector param json file
   -o pixcube     output h5 of pixcube
   -m             output sequential png files for movie
   --persistence  considering persistence (not supported yet!).
   --psf psf.fits  psf file (if not given, an analytic donuts model will be used.
""" 

from docopt import docopt             # command line interface
import numpy as np
from jis.pixsim import readflat as rf
from jis.pixsim import simpix_stable as sp
from jis.pixsim import gentraj
from jis.pixsim import makeflat as mf
from jis.jisplot import plotace 
from jis.photonsim.extract_json import mkDet
import tqdm
import astropy.io.fits as fits
import os
import sys
import h5py
import time
import matplotlib.pyplot as plt

# Constants. ##########
spixdim  = [32, 32] # subpixel dimension in a pixel (setting for intrapix pattern).


# Command line interface
if __name__ == '__main__':
    ts = time.time() # Starting time.

    args = docopt(__doc__)
    
    # Get parameters from command line
    tframe   = float(args['-s']) # Exposure time per frame (sec).
    nframe   = int(args['-f'])   # Number of frames.
    vd       = float(args['-d']) # Drift rate.
    det_json = args['--det']     # Det-param json filename.

    t = np.array(range(0,nframe))*tframe # Total time to simulate (sec).

    det = mkDet(det_json, spixdim=spixdim) # Making a detector class object.

    # Loading PSF
    if args['--psf']:
        phdul = fits.open(args['--psf'])
        psf = phdul[0].data
        psfheader = phdul[0].header
    else:
        psf = None
        
    # Loading light curve data.
    if args['-l']:
        lhdul = fits.open(args['-l'])
        ldata = lhdul[0].data
        lhead = lhdul[0].header
        tl    = ldata[0,:]
        injlc = ldata[1,:]
        if np.sum((t-tl)**2) > 0:
            sys.exit("time mismatch of "+str(args["-l"]))
    
    x_scale = float(args['-v'])   # stddev of ACE-x (ex. arcsec).
    y_scale = float(args['-w'])   # stddev of ACE-y (ex. arcsec).
    pix_scale = float(args['-p']) # pixel scale (ex. arcsec/pix).
    N = int(args['-n'])           # Number of pixels on a side in the output image.

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

    Tace = xhead['ACE-TOTT'] # Total time of the ACE data.
    Nace = len(xdata)

    #-----------------------------------------#

    # Setting artificial linear drift.
    if vd>0.0:
        drift_length  = vd*Tace
        #drift_azimuth=np.random.random()*2.0*np.pi
        drift_azimuth = np.pi/4.0
        drift_theta   = gentraj.gentraj_drift(Nace,drift_length,drift_azimuth)

    # Setting and plotting full trajectory.
    if vd>0.0:
        theta_full = np.array([xdata*x_scale/pix_scale+drift_theta[0,:], ydata*y_scale/pix_scale+drift_theta[1,:]])
    else:
        theta_full = np.array([xdata*x_scale/pix_scale, ydata*y_scale/pix_scale])
    ##### drift_theta should be divided by pix_scale (TK)?? #####

    plotace.trajectory(theta_full[0,:], theta_full[1,:])
    #    sys.exit()

    Nts_per_frame = int(tframe*Nace/Tace) # Number of timesteps per a frame.

    Nmargin  = 10
    Npixcube = int((np.max(np.abs(theta_full))+Nmargin)*2)
    pixdim   = [Npixcube, Npixcube] # adaptive pixel dimension in the aperture.

    # Preparing flat (intrapix/interpix). ##########
    ## Setting intrapixel pattern.
    intrapix = det.intrapix

    ## Setting GLOBAL interpixel pattern (flat).
    flat    = det.interpix
    gpixdim = np.shape(flat) # dimension for global pixel positions

    # Setting initial global position. #############
    x0 = (0.5*(np.shape(flat)[0]-spixdim[0]))
    y0 = (0.5*(np.shape(flat)[1]-spixdim[1]))
    x  = x0
    y  = y0

    # Setting persistence parameters. ##############
    if args["--persistence"]:
        tau = det.persistence['tau']
        rho = det.persistence['rho']
        # tau: detrapping timescales in sec.
        # rho: trapping fractions.

        NQc    = len(tau) #number of the defected pixel time scale
        Qtrap  = np.zeros((gpixdim[0],gpixdim[1],NQc)) #trapped charge
        xtau   = tframe/tau
        Qtsave = []

    # Making a LOCAL flat as interpix ##############
    jx, jy = np.int(x), np.int(y)
    interpix = rf.flat_interpix(flat, jx, jy, pixdim, figsw=0)

    lc = []        
    pixcube = np.zeros((Npixcube, Npixcube, nframe))
    for iframe in tqdm.tqdm(range(0,nframe)):

        #the positioning system is designed to match to the future upgrade of the large drift.
        #jx,jy=np.int(x),np.int(y)
        #interpix=rf.flat_interpix(flat,jx,jy,pixdim,figsw=0)
        
        # picking temporary trajectory and local position update
        istart = iframe*Nts_per_frame
        iend   = (iframe+1)*Nts_per_frame
        if iend >= Nace:
            print("insufficient time length of ACE fits.")
            sys.exit(-1)            

        # PSF center for iframe
        theta = np.copy(theta_full[:,istart:iend])
        theta = theta+np.array([pixdim]).T/2

        # Perform the PSF integration
        # output: array of images taken at each frame.
        pixar = sp.simpix(theta, interpix, intrapix, psf=psf)
        
        if args["--persistence"]:
            #persistence
            Qij = read_trapped_charge(Qtrap, jx, jy, pixdim, figsw=0)
            E0i = np.nansum(pixar,axis=2)
            Qij, Qi, Ei = persistence_const_array(xtau, rho, E0i, Qij)
            Qtrap = push_trapped_charge(Qtrap, jx, jy, pixdim, Qij)
            Qtsave.append(np.sum(Qtrap, axis=(0,1)))
            lctmp = np.mean(np.sum(Ei))
            lc.append(lctmp)

        pixcube[:,:,iframe] = np.sum(pixar,axis=2) # Making one exposure frame and storing it in pixcube.
            
    if args["-l"]:
        pixcube = pixcube*injlc

    te = time.time()   # End time.
    print(te-ts,"sec") # Show elapsed time.

    # Making output h5 data.
    with h5py.File(args["-o"], "w") as f:
        f.create_group("header")
        f.create_group("data")
        f.create_dataset("header/tframe", data=args["-s"])
        f.create_dataset("data/pixcube", data=pixcube)
        f.create_dataset("data/interpix", data=interpix)
    
    # Plotting light curve data.
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.plot(np.sum(pixcube[:,:,:],axis=(0,1))) # Light curve affected by interpix pattern.
    ax.plot(np.sum(pixcube.transpose()/interpix.transpose(),axis=(1,2))) # Light curve corrected for interpix pattern.
    ax.set_aspect(0.7/ax.get_data_ratio())
    plt.savefig("lc.png")

    # Saving each frame in pixcube.
    if args["-m"]:
        for iframe in range(0,nframe):
            plt.imshow(pixcube[:,:,iframe])
            plt.savefig("pixcube"+str(iframe)+".png")

