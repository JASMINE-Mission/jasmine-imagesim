import math
import numpy as np
import pyfftw
from scipy.ndimage import shift

def calc_psf(wfe, wN, k, WL, NP, Ntot, Stel, adata, M, aN, fN=520):
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
        k     (int)    : Number of wavelength data points.
        WL    (ndarray): Wavelength data (um).
        NP    (ndarray): Detected photon (electron) flux (e-/s/m^2/um).
        Ntot  (float)  : Total detected photon (electron) rate (e-/s/m^2).
        Stel  (float)  : Total area of the telescope aperture (m^2).
        adata (ndarray): Aperture mask data.
        M              : Inverse number of the PSF cell scale (mm/um).
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
        k, WL, NP, Ntot = extract_json.extsp(sp)

        # Reading aperture pattern.
        Stel, adata, aN, ahdr = readfits.read_aperture_mask("aperture.fits")

        # Reading wavefront error pattern.
        wfe,wN,whdr=readfits.read_wfe_map("wfe.fits")

        # Reading control parameters.
        cp = json.load(open("ctrl_param.json"))
        M = cp['M']['val']

        # Getting psf pattern.
        image = psf.calc_psf(wfe,wN,k,WL,NP,Ntot,Stel,adata,M,aN)
   

    """

    # Checking parities of wN and aN.
    # If their parities are different, the centers of wfe and apt will be different.
    if (wN+aN)%2 == 1:
        print("wN and aN should have the same parity")
        print("wN: {}; aN: {}".format(wN, aN))
        exit()
        
    wfer = np.nan_to_num( wfe )
    wfec = np.empty( (wN,wN),dtype='complex128')

    # 波長を k 個の点で規定している。
    # 口径 D で、計算領域の大きさを N とするとき、
    # フーリエ変換で得られる PSF は D/N x lambda/D=lambda/N の角度スケールになる。 
    # 異なるいくつかの波長 WL0, WL1, ..., WLn で生成したPSFのセルスケールを
    # 合わせようとおもったら、FFT計算領域を波長に比例させて
    # N=M WL (ここでWLはミクロン単位の波長)として設定し、
    # 得られた結果を PSF として fN x fN のサイズで切り取って加算すればよい。
    # なお、Mは計算する波長WLに対してNが整数になるように選んでおく。

    image = np.zeros( (fN, fN) ,dtype ='float' )

    # WL として、1.1, 1.2, ,,, 1.6 としているとき、
    # 1.1-1.2 の範囲のフォトンを 波長 1.15 で代表させて加え、
    # 1.2-1.3, ,,, 1.5-1.6 を加える、としていこう。 
    for i in range(k-1):
        WLm = (WL[i] + WL[i+1])/2 
        N = int(WLm*M)

        # initialize data for FFT
        data = pyfftw.zeros_aligned((N,N),dtype='complex128')

        #  put the mask on the data assuming N > aN
        i1 = int(N/2-aN/2)
        i2 = int(N/2+aN/2)
        data.real[i1:i2,i1:i2] = adata[:,:]
        
        # ここで、 WFE map を波長に反比例させて大きさをかえてから、
        # 位相成分として複素数化して data に掛ける
        wfec.real = np.cos(wfer/WLm)
        wfec.imag = np.sin(wfer/WLm)
        i3 = int(N/2-wN/2)
        i4 = int(N/2+wN/2)
        data[i3:i4,i3:i4] = data[i3:i4,i3:i4] * wfec
    
        # FFT
        ft  = pyfftw.interfaces.numpy_fft.fft2(data)
        fts = np.fft.fftshift(ft)

        # 結果を入射光子数(電子数)の重みを付けて加算する。
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
    
    # 最後に規格化しておく、値は1秒あたりのelectron数になるように
    s     = np.sum(image)
    image = image/s * Ntot * Stel

    return image
