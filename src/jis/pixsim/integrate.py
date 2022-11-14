import numpy as np
from jis.pixsim.addnoise import mk_readnoise
from jis.pixsim.addnoise import mk_shotnoise
import matplotlib.pylab as plt

def integrate(pixar, jx, jy, texp, dt, det, raw=False, addnoise=True, digitize=True):
    """
    Summary:
        This function integrates the pixar array
        calculated by simpix with taking the sampling
        timing. This function also adds readnoise and
        shotnoise, and digitize. Saturation is also considered.

    Args:
        pixar  (ndarray): Movie data created by simpix.
                          The shape is (X, Y, Z) not (Z, Y, X).
        jx, jy (int)    : Global pixel positions correspond to the
                          origin of the simulated local image (pixar).
        texp   (float)  : Exposure time in sec.
        dt     (float)  : Timestep of the simulated data (pixar).
        det    (detctor): Detector class object.
        raw    (bool)   : Provide two images before subtraction?
        addnoise (bool) : Switch of noise-addition function.
        digitize (bool) : Switch of digitization process.

    Returns:
        adu (ndarray): Integrated data.
                       If raw=True, [adu2, adu1] are given.
                       If raw=False, the output is adu2-adu1.
                       If digitiza=True, the unit is adu.
                       If digitize=False, the unit is e-.
    
    """

    pixar = pixar.T
    nz, ny, nx = pixar.shape

    time = (np.indices(pixar.shape)[0]+1) * dt

    # Local coordinate.
    y_local, x_local = np.indices([ny, nx])

    # Global coordinate.
    y_global = y_local + jy
    x_global = x_local + jx

    # Channel coordinate.
    y_ch = np.copy(y_global)
    x_ch = x_global%det.readparams.ncol_ch

    # First sampling time since the end of reset.
    t1 = det.readparams.t_overhead + \
         det.readparams.tsmpl*(det.readparams.npix_read_per_row*y_ch+det.readparams.npix_pre+x_ch)

    # Second sampling time since the end of reset.
    t2 = t1 + texp

    # Weight array for the first sampling.
    tmp = (time-t1)/det.readparams.tsmpl
    w1 = np.zeros(shape=pixar.shape)
    w1[tmp<=0.] = 1.0
    w1[tmp>=1.] = 0.0
    pos = np.where((tmp>0.)*(tmp<1.))
    w1[pos] = 1.0 - tmp[pos]

    # Weight array for the second sampling.
    tmp = (time-t2)/det.readparams.tsmpl
    w2 = np.zeros(shape=pixar.shape)
    w2[tmp<=0.] = 1.0
    w2[tmp>=1.] = 0.0
    pos = np.where((tmp>0.)*(tmp<1.))
    w2[pos] = 1.0 - tmp[pos]

    # Integration (ideal).
    integ1 = np.sum(pixar*w1, axis=0)
    integ2 = np.sum(pixar*w2, axis=0)

    # Adding shotnoise.
    shotnoise1, seed1 = mk_shotnoise(integ1)
    shotnoise12, seed2 = mk_shotnoise(integ2-integ1)

    signal1 = integ1 + shotnoise1
    signal2 = signal1 + (integ2 - integ1) + shotnoise12

    ###################################################
    # If needed, non-linearity can be implemented here
    # instead of the saturation cut below.
    ###################################################

    # Saturation cut
    signal1[signal1>det.fullwell] = det.fullwell
    signal2[signal2>det.fullwell] = det.fullwell

    # Adding readnoise.
    if addnoise:
        adu1 = signal1 + mk_readnoise(signal1.shape, det.readnoise)[0]
        adu2 = signal2 + mk_readnoise(signal2.shape, det.readnoise)[0]
    else:
        adu1 = signal1
        adu2 = signal2
    # Digitization.
    if digitize:
        adu1 = np.round(adu1/det.gain)
        adu2 = np.round(adu2/det.gain)

    if raw:
        return adu2.T, adu1.T
    else:
        # Subtraction.
        adu = (adu2 - adu1).T
        return adu
