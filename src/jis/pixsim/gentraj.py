"""Generating the trajectory array

"""

import numpy as np
import argparse
import pylab
import matplotlib.pyplot as plt


def gentraj_drift(Nts, drift_length, drift_azimuth):
    """generate a linear pointing trajectory in detpix

    Args:
        Nts (int): number of the trajectory array 
        drift_length (float): length of the linear trajectory in detpix
        drift_azimuth (float): azimuth angle (radian) of the linear trajectory  

    Returns:
        nd array: trajectory array in detpix
    """
    thetax = np.linspace(0, drift_length, Nts) * np.cos(drift_azimuth)
    thetay = np.linspace(0, drift_length, Nts) * np.sin(drift_azimuth)
    theta = np.array([thetax, thetay])
    return theta


def gentraj_random(ntime, basepos, nsub, basesig=1.0, subsig=0.1, seed=None):
    """Random Gaussian trajectory generator 

    Args:
       ntime: the number of total time bin
       basepos: base position
       nsub: number of the time bin at  one orbit, 
       basesig: sigma for the initial pointing accuracy
       subsig: sigma during one orbit accuracy

    """
    ndata = np.int(ntime / nsub)
    np.random.seed(seed)
    thetax = np.array([])
    thetay = np.array([])
    retposx = np.random.randn(ndata) * basesig
    retposy = np.random.randn(ndata) * basesig
    for i in range(0, ndata):
        thetax = np.hstack(
            [thetax, basepos[0] + retposx[i] + np.random.randn(nsub) * subsig])
        thetay = np.hstack(
            [thetay, basepos[1] + retposy[i] + np.random.randn(nsub) * subsig])

    theta = np.array([thetax[0:ntime], thetay[0:ntime]])
    np.random.seed()

    return theta


if __name__ == "__main__":

    nsub = 16
    ntime = 1024 * nsub
    basepos = [6, 6]
    pixsec = 29.0  #arcsec/pix
    sigsec_interframe = 0.0
    sigsec_intraframe = 11.0  #arcsec

    theta = gentraj_random(ntime,
                           basepos,
                           nsub,
                           basesig=sigsec_interframe / pixsec,
                           subsig=sigsec_intraframe / pixsec,
                           seed=None)
    fig = plt.figure(figsize=(10, 10))
    plt.gca().invert_yaxis()
    plt.plot(theta[0, ::], theta[1, :], ".", c="green", markersize=2)
    plt.plot(theta[0, ::nsub], theta[1, ::nsub], "+", c="red", lw=1)
    plt.show()
