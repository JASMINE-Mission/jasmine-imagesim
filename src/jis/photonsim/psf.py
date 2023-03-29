import math
import numpy as np
import pyfftw
from scipy.ndimage import shift

def calc_psf(wfe, wN, WL, NP, Ntot, Stel, adata, M, aN, fN=520):
    """
    This function calculates the psf in e-/sec/fp-cell
    based on the wavefront error (wfe),
    the spectral information (NP), and
    the aperture mask data (adata).
    The fp-cell scale is defined by the parameter M as
    (1/M) x 10^-3 rad/fp-cell.

    Args:
        wfe   (ndarray): Wavefront error data (um).
        wN    (int)    : Number of fp-cells on a side of the wfe data.
        WL    (ndarray): Wavelength data (um).
        NP    (ndarray): Detected photon (electron) flux (e-/s/m^2/um).
        Ntot  (float)  : Total detected photon (electron) rate (e-/s/m^2).
        Stel  (float)  : Total area of the telescope aperture (m^2).
        adata (ndarray): Aperture mask data.
        M     (int)    : Inverse number of the PSF cell scale (mm/um).
        aN    (int)    : Number of apt-cells of the aperture mask data.
        fN    (int)    : Number of fp-cells of the output psf data (Default: 520).

    Returns:
        image (nadarray): fN x fN fp-cell array of the psf (e-/s/fp-cell).
                          The fp-cell scale is (1/M) x 10^-3 rad/fp-cell.

    Example:
        import json
        from jis.photonsim import extract_json
        from jis.photonsim import readfits

        # Reading source information.
        sp = json.load(open("source.json"))
        WL, NP, Ntot = extract_json.extSp(sp)

        # Reading aperture pattern.
        Stel, adata, aN, ahdr = readfits.read_aperture_mask("aperture.fits")

        # Reading wavefront error pattern.
        wfe,wN,whdr=readfits.read_wfe_map("wfe.fits")

        # Reading control parameters.
        cp = json.load(open("ctrl_param.json"))
        M = cp['M']['val']

        # Getting psf pattern.
        image = psf.calc_psf(wfe,wN,WL,NP,Ntot,Stel,adata,M,aN)


    """

    # Checking parities of wN and aN.
    # If their parities are different, the centers of wfe and apt will be different.
    if (wN+aN)%2 == 1:
        print("wN and aN should have the same parity")
        print("wN: {}; aN: {}".format(wN, aN))
        exit()

    wfer = np.nan_to_num( wfe )
    wfec = np.empty( (wN,wN),dtype='complex128')

    # We define the wavelength grid with k points.
    # When the aperture diameter is D, and the size of the calculation region is N,
    # the PSF image obtained through FFT will have an angle grid scale of
    # D/N x lambda/D=lambda/N. To equalize the PSF angle scales for different
    # wavelengths, WL0, WL1, ..., WLn, it is a good way to vary the calculation
    # region N in prop. to wavelength as, N=M WL (WL: wavelength in um) and
    # add the obtained PSF images after extracting the central fN x fN region.
    # For accuracy, M should be chosen such that N is an integer for all wavelengths.

    image = np.zeros( (fN, fN) ,dtype ='float' )

    # We calculate the photons in a certain range between two wavelengths
    # as those at the mean wavelength of the two and add them togather.
    for i in range(len(WL)-1):
        WLm = (WL[i] + WL[i+1])/2
        N = int(WLm*M)

        # initialize data for FFT
        data = pyfftw.zeros_aligned((N,N),dtype='complex128')

        #  put the mask on the data assuming N > aN
        i1 = int(N/2-aN/2)
        i2 = int(N/2+aN/2)
        data.real[i1:i2,i1:i2] = adata[:,:]

        # Converting the wfe data to phase error and
        # multiply the complex data created from the phase error
        # to the aperture pattern.
        wfec.real = np.cos(wfer/WLm*2.*np.pi)
        wfec.imag = np.sin(wfer/WLm*2.*np.pi)
        i3 = int(N/2-wN/2)
        i4 = int(N/2+wN/2)
        data[i3:i4,i3:i4] = data[i3:i4,i3:i4] * wfec

        # FFT
        ft  = pyfftw.interfaces.numpy_fft.fft2(data)
        fts = np.fft.fftshift(ft)

        # Weighting with the number of photons/electrons and
        # adding together.
        i1 = int(N/2-fN/2)
        i2 = int(N/2+fN/2)
        NPm   = (NP[i] + NP[i+1])/2

        tmp_img = NPm*(fts.real[i1:i2, i1:i2]**2.+fts.imag[i1:i2, i1:i2]**2.)
        offset = 0.0
        if   (N%2==0) and (fN%2==1):
            offset=1.0
        elif (N%2==0) and (fN%2==0):
            offset=0.5
        elif (N%2==1) and (fN%2==0):
            offset=0.5
        tmp_img = shift(tmp_img, [-offset, -offset], order=1)

        image = image + tmp_img

    # Normalizing the resutl such that the total value equals to the electron rate.
    s     = np.sum(image)
    image = image/s * Ntot * Stel

    return image


def calc_dummy_psf(wfe, wN, WL, NP, Ntot, Stel, adata, M, aN, fN=520):
    """
    This function calculates a dummy point-spread function image in e-/sec/fp-cell.
    The arguments are the same as `calc_psf` but not used except for `Ntot`, `Stel`, and `fN`.
    Currently the FWHM of the PSF is the one-tenth of the fp-array.

    Args:
        wfe   (ndarray): Wavefront error data (um).
        wN    (int)    : Number of fp-cells on a side of the wfe data.
        WL    (ndarray): Wavelength data (um).
        NP    (ndarray): Detected photon (electron) flux (e-/s/m^2/um).
        Ntot  (float)  : Total detected photon (electron) rate (e-/s/m^2).
        Stel  (float)  : Total area of the telescope aperture (m^2).
        adata (ndarray): Aperture mask data.
        M     (float)  : Inverse number of the PSF cell scale (mm/um).
        aN    (int)    : Number of apt-cells of the aperture mask data.
        fN    (int)    : Number of fp-cells of the output psf data (Default: 520).

    Returns:
        image (nadarray): fN x fN fp-cell array of the psf (e-/s/fp-cell).
                          The fp-cell scale is (1/M) x 10^-3 rad/fp-cell.
    """
    fwhm = 2.0
    sigma = fwhm/2.0/np.sqrt(2*np.log(2))
    arr = np.linspace(-10,10,fN)
    xx,yy = np.meshgrid(arr,arr)
    image = np.exp(-(xx**2+yy**2)/2.0/sigma**2)
    return image/image.sum() * Ntot * Stel


def calc_gauss_psf(fwhm, Ntot, Stel, M, fN=520):
    """
    This function creates a gaussian PSF image in e-/sec/fp-cell.

    Args:
        fwhm (float): PSF FWHM (rad).
        Ntot (float): Total detected photon (electron) rate (e-/s/m^2).
        Stel (float): Total area of the telescope aperture (m^2).
        M    (float): Inverse number of the PSF cell scale (mm/um).
        fN   (int)  : Number of fp-cells of the output psf data (Default: 520).

    Returns:
        image (ndarray): fN x fN fp-cell array of the psf (e-/s/fp-cell).
                         The fp-cell scale is (1/M) x 10^-3 rad/fp-cell.
    """
    fp_cell_scale = 1.e-3/M # rad/fp-cell.

    fwhm_fp_cell  = fwhm/fp_cell_scale
    sigma         = fwhm_fp_cell/2./np.sqrt(2.*np.log(2.))

    arr = np.arange(0,fN)-(fN-1.)/2.
    xx, yy = np.meshgrid(arr, arr)
    image  = np.exp(-(xx**2.+yy**2.)/2./sigma/sigma)
    return image/image.sum()*Ntot*Stel
