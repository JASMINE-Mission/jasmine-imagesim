import numpy as np
import pylab 
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import argparse
import shutil
import readflat as rf
import simpixlc as sp
import gentraj as gt
import addnoise as an
import transitmodel as tm
import pld

import sys

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-s', nargs=1, default=[0.5],help='sigma of PSF gaussian. typical Kepler = 0.5 ', type=float)
    parser.add_argument('-b', nargs=1, default=[1800.0], help='second per bin', type=float)
    parser.add_argument('-t', nargs=1, default=["eck2like"],help='choose trajectory: k2like, random, ranwalk', type=str)
    parser.add_argument('-p', nargs=2, default=[8.5,7.5], help='base pix position', type=float)
    parser.add_argument('-g', nargs=1, default=[2.0], help='sigma trajectory (basesig or dr or radius)', type=float)
    parser.add_argument('-q', nargs=1, default=[500.0], help='photon noise (ppm/bin)', type=float)
    parser.add_argument('-d', nargs='+', default=[0.001,0.002,0.003], help='inject transit depth', type=float)
    parser.add_argument('-z', nargs='+', default=[436,636,836], help='inject transit center', type=int)
    parser.add_argument('-i', nargs=1, default=[0], help='tag', type=int)
    parser.add_argument('-r', nargs=1, default=[None], help='random seed', type=int)
    parser.add_argument('-n', nargs=1, default=[4022], help='# of targ images (ndata)', type=int)
    parser.add_argument('-x', nargs=1, default=[180], help='subsampling factor. # of frames/readout', type=int)
    #parser.add_argument('-l', nargs='+', default=[1e0, 1e5, 1e10, 1e15], help='regularization parameter', type=float)
    parser.add_argument('-l', nargs='+', default=[1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14], help='regularization parameter', type=float)
    det=True
    ref=True

    args = parser.parse_args()
    sigma= args.s[0]
    sigma2=sigma*sigma
    bin2sec=args.b[0]
    trajtype=args.t[0]
    ppx,ppy=args.p
    bsig=args.g[0]
    ppm=args.q[0]
    tdepth=args.d
    t0=args.z
    seed=args.r[0]
    ndata = args.n[0]
    nsub = args.x[0]
    ntime=ndata*nsub
    ntransit=len(tdepth)
    itag = '{0:04d}'.format(args.i[0])

    print("total frame:",ntime)
    print("# of data:", ndata)

    ### FILE I/O ###
    dirh="/home/kawahara/virtual-jasmine/pixsim/data/"
    outdir="out"+str(itag)+"."+trajtype+"_bsig"+str(bsig)+"_psf"+str(sigma)+"_ppm"+str(np.int(ppm))

    filetemp="template.fits.gz"
    filet="pixel.fits.gz"
    outpath=os.path.join(dirh,outdir)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    shutil.copy(os.path.join(dirh,filetemp),os.path.join(outpath,filet))

    ### PIXEL Settings ###
    #pixdim = [25, 25]  # pixel dimension in the aperture
    pixdim=[18,16] # pixel dimension in the aperture
    spixdim=[32,32] # subpixel dimension in a pixel

    
    #trajectory
    if trajtype == "k2like":
        theta=gt.gentraj_k2like(ntime,np.array([ppx,ppy]),nsub,bsig,seed=seed)
    elif trajtype == "random":
        theta = gt.gentraj_random(ntime, np.array([ppx,ppy]),nsub, bsig,seed=seed)
    elif trajtype == "ranwalk":
        theta = gt.gentraj_ranwalk(ntime, np.array([ppx,ppy]), bsig,seed=seed)
    elif trajtype == "ecrandom":
        theta = gt.gentraj_encircled_random(ntime, np.array([ppx,ppy]), radius=bsig,seed=seed)
    elif trajtype == "eck2like":
        theta = gt.gentraj_encircled_k2like(ntime, np.array([ppx,ppy]),nsub, radius=bsig, lpix=1.4, gradsig=3.0,seed=seed)
    else:
        print("No trajectory type. Use k2like.")
        sys.exit()

    #interpixel
    np.random.seed(seed)
    jx,jy=np.random.rand(2)
    #jx=0.5
    #jy=0.5
    hayapix=1024
    jx=np.int(jx*(1024-spixdim[0]))
    jy=np.int(jy*(1024-spixdim[1]))
    interpix=rf.haya2interpix(jx,jy,pixdim)

    #intrapixel
    filex="intravx.csv"
    filey="intravy.csv"
    intrapix=rf.read_intrapix(filex,filey,spixdim)

    ### COPY ###
    #hdulist = fits.open(os.path.join(outpath, filet))
    #cp0=hdulist[0].header
    #cp1h=hdulist[1].header
    #cp1d=hdulist[1].data
    #hdulist.close()
    ### PASTE ###
    #print(np.shape(cp1d))
    #print(np.shape(cp1d["TIME"]))
    #fileo="test.fits"
    #hdulist = fits.open(os.path.join(outpath, fileo))
    #primary=fits.PrimaryHDU([])
    #tt=fits.BinTableHDU([])
    #hdulist=fits.HDUList([primary,tt])
    #hdulist[0].header=cp0
    #hdulist.writeto(os.path.join(outpath, fileo),overwrite=True)
    #sys.exit()

    #generate
    pixar=sp.simpix(theta,interpix,intrapix,sigma2=sigma2)
    print(np.shape(pixar))
    pixar=np.sum(pixar.reshape(pixdim[0],pixdim[1],ndata,nsub),axis=3)
    print(np.shape(pixar))
    #sys.exit()

    cfac=an.convfactor_photon(ppm,np.mean(np.sum(pixar,axis=(0,1))))
    for itransit in range(0,ntransit):
        pixar=pixar*tm.inject_transit(ndata, t0[itransit], tdepth[itransit], ndata)
    pixar = np.random.poisson(cfac*pixar)

    #reference (almost no systematic)
    if ref==True:
        interpixref = np.ones(pixdim)
        intrapixref = np.ones(spixdim)
        pixarref = sp.simpix(theta, interpixref, intrapixref, sigma2=sigma2)
        pixarref=np.sum(pixarref.reshape(pixdim[0],pixdim[1],ndata,nsub),axis=3)
        cfacref = an.convfactor_photon(ppm, np.mean(np.sum(pixarref, axis=(0, 1))))
        for itransit in range(0, ntransit):
            pixarref = pixarref * tm.inject_transit(ndata, t0[itransit], tdepth[itransit], ndata)
        pixarref = np.random.poisson(cfacref * pixarref)

    ### UPDATING ###
    hdulist = fits.open(os.path.join(outpath,filet),mode="update")
    testd=hdulist[1].data
    testd["RAW_CNTS"]=np.copy(pixar.T)
    testd["FLUX"]=np.copy(pixar.T/(bin2sec))
    testd["FLUX_ERR"]=np.copy(np.sqrt(pixar.T/(bin2sec)))

    ### TDIM4,9 ###
    #testh = hdulist[1].header
    #for tag in ["TDIM4","TDIM5","TDIM6","TDIM7","TDIM8","TDIM9"]:
    #    inpt="("+str(pixdim[0])+","+str(pixdim[1])+")"
    #    testh[tag]=inpt
    hdulist.close()

    ### COPY ###
    hdulist = fits.open(os.path.join(outpath, filet))
    cp0=hdulist[0]
    cp1h=hdulist[1].header
    cp1d=hdulist[1].data
    hdulist.close()
    ### PASTE ###

    if det==True:
        ### DETREND BY EVEREST(MOD) ###
        detflux=pld.detrend(os.path.join(outpath, filet),lamin=args.l)


    ### SAVE INFO ###
    np.savetxt(os.path.join(dirh,outdir,"pixelpos.txt"),theta.T)
    info="Trajectory type="+str(trajtype)+"\n"\
    +"PSF sigma="+str(sigma)+"pix\n"\
    +"1bin="+str(bin2sec)+"sec\n"\
    +"base position X,Y="+str(ppx)+","+str(ppy)+"\n"\
    +"base sigma="+str(bsig)+"\n"\
    +"photon noise/bin="+str(ppm)+"ppm\n"\
    +"flat position="+str(jx)+","+str(jy)+"\n"
    f=open(os.path.join(dirh,outdir,"log.txt"),"w")
    f.write(info)
    f.close()


    ### FIGURES ###
    fig = plt.figure(figsize=(13, 3))
    i = np.int(ndata / 2)
    ax = fig.add_subplot(121)
    a = plt.imshow(pixar[:, :, i].T, interpolation="nearest")
    plt.colorbar(a)
    plt.gca().invert_yaxis()
    ax2 = fig.add_subplot(122)
    plt.plot(np.sum(pixar, axis=(0, 1)) / np.mean(np.sum(pixar, axis=(0, 1))) - 1, ".")
    plt.ylabel("integrated LC")
    plt.title("std=" + str(np.std(np.sum(pixar, axis=(0, 1)) / np.mean(np.sum(pixar, axis=(0, 1))) - 1)))
    plt.xlabel("time step")
    plt.savefig(os.path.join(dirh, outdir, "example.png"))

    fig = plt.figure(figsize=(4, 4))
    plt.imshow(interpix.T)
    plt.gca().invert_yaxis()
    plt.plot(theta[0, :], theta[1, :],".", c="white", lw=1)
    plt.savefig(os.path.join(dirh, outdir, "flattraj.png"))

#    np.savez(os.path.join(dirh, outdir, "pixar"),pixar)
    
    ### DETREND
    if det==True:

        fig = plt.figure(figsize=(12,6))
        ax = fig.add_subplot(311)
        lc=np.sum(pixar, axis=(0, 1))/bin2sec
        ind = np.array(range(0, len(lc)))
        mask=(ind > -1)
        for itransit in range(0,ntransit):
            mask = mask*(np.abs(ind - t0[itransit]) > 30)

        ax.plot(lc, ".")
        plt.ylabel("Raw Data (e-/s)")
        ppmraw = np.std(lc[mask]) / np.mean(lc[mask]) * 1.e6
        ax.annotate(str(ppmraw)+" ppm/bin", \
        xy=(0, 0.8), xycoords='axes fraction', fontsize=12)
        np.savez(os.path.join(dirh, outdir, "rawlc.npz"),lc)

        ax2 = fig.add_subplot(312)
        lc=detflux
        ax2.plot(lc, ".")
        plt.ylabel("Detrended (e-/s)")
        ppmdet = np.std(lc[mask]) / np.mean(lc[mask]) * 1.e6
        ax2.annotate(str(ppmdet)+" ppm/bin", \
        xy=(0, 0.8), xycoords='axes fraction', fontsize=12)
        np.savez(os.path.join(dirh, outdir, "detlc.npz"), lc)

        ax3 = fig.add_subplot(313)
        lc=np.sum(pixarref, axis=(0, 1))/bin2sec
        ax3.plot(lc, ".")
        ax3.plot(ind[mask],lc[mask], ".",color="green")
        plt.ylabel("No fluctuation (e-/s)")
        ppmref=np.std(lc[mask])/np.mean(lc[mask])*1.e6
        ax3.annotate(str(ppmref)+" ppm/bin", \
        xy=(0, 0.8), xycoords='axes fraction', fontsize=12)
        np.savez(os.path.join(dirh, outdir, "ref.npz"), lc)
        plt.savefig(os.path.join(dirh, outdir, "comparison.png"))

        ### SAVE INFO ###
        np.savetxt(os.path.join(dirh, outdir, "detrend_result.txt"),\
                   np.array([ppm,ppmraw,ppmdet,ppmref]))
        #plt.show()
