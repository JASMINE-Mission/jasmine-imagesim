import pytest
import numpy as np
from jis.pixsim.simpix_stable import simpix
from jis.photonsim.extract_json import Detector
from jis.pixsim import readflat as rf

def test_simpix(fig=False):
    filename_detjson="det_fortest.json"
    spixdim  = [16,16] # subpixel dimension in a pixel (setting for intrapix pattern).
    pixdim   = [16,16] # adaptive pixel dimension in the aperture.
    det = Detector.from_json(filename_detjson)
    intrapix = det.flat.intrapix
    flat    = det.flat.interpix
    gpixdim = np.shape(flat) # dimension for global pixel positions
    print(flat)
    # Setting initial global position. #############
    x0 = (0.5*(np.shape(flat)[0]-spixdim[0]))
    y0 = (0.5*(np.shape(flat)[1]-spixdim[1]))
    x  = x0
    y  = y0
    jx, jy = np.int32(x), np.int32(y)
    interpix = rf.flat_interpix(flat, jx, jy, pixdim, figsw=0)
    theta=np.array([[6.5],[6.5]])
    
    val=simpix(theta, interpix, intrapix)[:,:,0]
    ##np.savez("ref.npz",val)
    if fig:
        import matplotlib.pyplot as plt
        c=plt.imshow(val)
        plt.colorbar(c)
        plt.show()
        
    valref=np.load("ref.npz")["arr_0"]
    assert np.sum((val-valref)**2)<1.e-16 


if __name__ == "__main__":
    print(test_simpix())
