import numpy as np
import pylab
import matplotlib.pyplot as plt
import os, pkgutil
from io import BytesIO
from astropy.io import fits
import argparse

def flat_interpix(flat, ix, iy, pixdim, figsw=0):
    """
    Summary:
        This function calculates an interpixel flat pattern
        based on the given flat data 'flat'. From the flat data,
        this function selects a region defined by ix, iy, and pixdim
        as x = ix to (ix+pixdim[0]-1) and y = iy to (iy+pixdim[1])
        (these pixel numbers are 0-origin values).
        After normalizing it with the mean value, this function
        returns the result as a 2-d array data.

        If figsw is set to 1, a preview is shown.

    Args:
        flat   (array): Original 2-d array data.
                        0th axis is x axis, 1st axis is yaxis.
        ix, iy (int)  : Initial pixel indices to select a region
                        from the flat data (0 origin).
        pixdim (array): Pixel dimensions (2-values).
                        The 1st/2nd values are the x/y dimensions.
        figsw  (int)  : Figure switch (Default: 0).
                        If figsw=1, a preview figure is shown.

    Returns:
        sflat  (array): Selected and normalized interpixel pattern.

    """

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


def haya2interpix(ix, iy, pixdim,
                  filename="pixsim/data/v_flatn.fits",figsw=0):
    """
    Summary:
        This function returns a interpixel flat pattern
        based on the haya2 data.
        (A wrapper of flat_interpix for the haya2 data).

    Args:
        ix, iy (int)  : Initial pixel indices to select a region
                        from the flat data (0 origin).
        pixdim (array): Pixel dimensions (2-values).
                        The 1st/2nd values are the x/y dimensions.
        filename (str): Name of the haya2 data (Default: pixsim/data/v_flatn.fits).
        fitsw (int)   : Figure switch (Default: 0).
                        If figsw=1, a preview figure is shown.

    Returns:
        sflat (array): Selected and normalized flat data.

    """
    try:
        data = pkgutil.get_data('jis', filename)
        hdulist = fits.open(BytesIO(data))
    except Exception as e:
        hdulist = fits.open(filename)

    flat=hdulist[0].data

    sflat = flat_interpix(flat, ix, iy, pixdim, figsw=figsw)

    return sflat

def read_intrapix(filex, filey, spixdim):
    """
    Summary:
        This function returns an intrapixel flat pattern.
        The pattern is read from files given as `filex` and `filey`,
        that should contain one-dimensional intrapix patterns along with
        x- and y-axis, respectively. The function first searches the files
        in the package directory, and then looks for them as relative paths.
        The flat pattern is defined as the product of the both flat patterns.
        The numbers of sub-pixels are set by `spixdim`.

    Args:
        filex   (str)  : Filename of x-direction intrapixel pattern.
        filey   (str)  : Filename of y-direction intrapixel pattern.
        spixdim (array): 2-value array of the numbers of sub-pixels.
                         (1st/2nd values are numbers in x/y directions)

    Returns:
        intrapix (array): Calculated and normalized interpixel pattern.
                          (Normalized by the mean value)
    """

    try:
        datax = pkgutil.get_data('jis', filex)
        datay = pkgutil.get_data('jis', filey)
        xd = np.loadtxt(BytesIO(datax), delimiter=',')
        yd = np.loadtxt(BytesIO(datay), delimiter=',')
    except Exception as e:
        xd = np.loadtxt(filex,delimiter=",")
        yd = np.loadtxt(filey,delimiter=",")

    pxvals=np.linspace(-0.5,0.5,spixdim[0])
    xintra=np.interp(pxvals, xd.T[0],xd.T[1])

    pyvals=np.linspace(-0.5,0.5,spixdim[1])
    yintra=np.interp(pyvals, yd.T[0],yd.T[1])

    intrapix=xintra*np.array([yintra]).T
    intrapix=intrapix/np.mean(intrapix)

    return intrapix
