import pytest
import numpy as np
from jis.pixsim.gentraj import gentraj_random

def test_gentraj_random():

    ref=np.array([[5.79965899,5.5930119,6.32825807,5.12700256,6.66182515,5.71126635,6.12101483,5.90541124,6.55459267,5.21856732,5.87770382,5.85432421,6.43005048,5.58279986,5.9345962,5.66701922],[6.01601211,6.22106784,5.58252376,6.43420555,6.34198269,6.1906013,6.34170398,5.74065495,5.95338647,5.64505297,5.89838728,6.20116931,5.73764592,5.84950728,5.73934829,5.67940476]])
    ntime=16
    basepos=[6, 6]
    pixsec=29.0 #arcsec/pix
    sigsec_interframe=0.0
    sigsec_intraframe=11.0 #arcsec

    theta=gentraj_random(ntime,basepos,ntime,basesig=sigsec_interframe/pixsec,subsig=sigsec_intraframe/pixsec,seed=1)
    assert np.sum((theta-ref)**2)<3.e-16

if __name__ == "__main__":
    print(test_gentraj_random())


