import numpy as np
import argparse
import readflat as rf
import os
import simpix_stable as sp

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This is just simple test code to generate light curves from the pixel-level simulation. You need to specify the trajectory in test.csv (or file called by -t option).')
    parser.add_argument('-t', nargs=1, default=["test.csv"],help='trajectory file', type=str)
    parser.add_argument('-s', nargs=1, default=[0.5],help='sigma of PSF gaussian. typical Kepler = 0.5 ', type=float)

    args = parser.parse_args()
    sigma= args.s[0]
    sigma2=sigma*sigma
    trajfile=args.t[0]

    pixdim=[32,32] # pixel dimension in the aperture
    spixdim=[32,32] # subpixel dimension in a pixel
    ntime=1024
    
    #intrapixel
    filex="intravx.csv"
    filey="intravy.csv"
    datapath="../data/intrapix"
    intrapix=rf.read_intrapix(filex,filey,spixdim,datapath)
    
    #interpixel from the Hayabusa2 ONC flat

    jx,jy=np.random.rand(2)
    hayapix=1024
    jx=np.int(jx*(1024-spixdim[0]))
    jy=np.int(jy*(1024-spixdim[1]))
    interpix=rf.haya2interpix(jx,jy,pixdim)

    #trajectory
    theta=np.loadtxt(trajfile).T
    print(np.shape(theta))
    #generate
    pixar=sp.simpix(theta,interpix,intrapix,sigma2=sigma2)

    ### FIGURES
    import matplotlib.pyplot as plt
    import pylab
    for i in range(0,np.shape(pixar)[2]):
        fig = plt.figure()
        plt.imshow(pixar[:,:,i])
        plt.show()
