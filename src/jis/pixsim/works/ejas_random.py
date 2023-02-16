import numpy as np
import argparse
import readflat as rf
import os
import simpix_stable as sp
import makeflat as mf
import matplotlib.pyplot as plt
import pylab
import transitmodel
import sys
#
# X,Y in input/output are reverted in y, x in kernel
#


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='This code generates the light curve with a planet assuming a random selection of a position in a detector plane. No persistence. Interpix fluctuation is 1 percent gaussian. Intrapix model from CCD.')
    parser.add_argument('-t', nargs=1, default=["test.csv"],help='trajectory file', type=str)
    parser.add_argument('-s', nargs=1, default=[0.5],help='sigma of PSF gaussian. typical Kepler = 0.5 ', type=float)
    parser.add_argument('-f', help='flat correction?', action='store_true')
    args = parser.parse_args()
    sigma= args.s[0]
    sigma2=sigma*sigma
    trajfile=args.t[0]

    pixdim=[20,20] # pixel dimension in the aperture
    spixdim=[32,32] # subpixel dimension in a pixel
    ntime=1024

    #generate Gaussian flat 
    flat = mf.gaussian_flat(sigma=0.01)
    #intrapixel
    filex="intravx.csv"
    filey="intravy.csv"
    dirh="../data/intrapix"
    intrapix=rf.read_intrapix(filex,filey,spixdim,dirh)

    #trajectory
    sigmaXY=0.567
    theta_cen = np.array(pixdim)/2.0
    ndisturb=256
    nframe = 3000
    t=np.array(range(0,nframe))*7.0 #sec

    #planet model
    Rpin=1.0
    injlc, b=transitmodel.gentransit(t/3600/24,t0=180.0/60/24,Rp=Rpin)

    lc=[]
    for iframe in range(0,nframe):

        #interpixel
        jx,jy=np.random.rand(2)
        jx=int(jx*(np.shape(flat)[0]-spixdim[0]))
        jy=int(jy*(np.shape(flat)[1]-spixdim[1]))
        interpix=rf.flat_interpix(flat,jx,jy,pixdim,figsw=0)

        theta=  np.random.normal(0.0,sigmaXY,2*ndisturb).reshape((ndisturb,2)) + np.array([theta_cen])
        theta = theta.T
        pixar=sp.simpix(theta,interpix,intrapix,sigma2=sigma2)
        
        if args.f:
            #flat correction==========
            est_error = 0.001 #estimation error of flat
            error=(np.random.normal(0.0,est_error,pixdim[0]*pixdim[1]).reshape(pixdim[0],pixdim[1]))*np.mean(interpix)
            interpix_est=interpix+error
            for i in range(0,np.shape(pixar)[2]):
                pixar[:,:,i] = pixar[:,:,i]/interpix_est
            #=======================================
        
        lctmp=np.mean(np.sum(pixar,axis=(0,1)))
        lc.append(lctmp)

    lc=lc/np.mean(lc)*injlc


    Nbin=43
    lcbin=[]
    tbin=[]
    for i in range(0,int(len(lc)/Nbin)):
        lcbin.append(np.mean(lc[i*Nbin:(i+1)*Nbin]))
        tbin.append(np.mean(t[i*Nbin:(i+1)*Nbin]))
    lcbin=np.array(lcbin)
    tbin=np.array(tbin)
        
    print("STD=",np.std(lcbin))
    plt.plot(t/60,lc,".",alpha=0.3,label="frame")
    plt.plot(tbin/60,lcbin,"s",alpha=1.0,label="5 min binning")
    plt.legend()
    plt.xlabel("[min]")
    if args.f:
        plt.savefig("lc.flat.png")
    else:
        plt.savefig("lc.png")
