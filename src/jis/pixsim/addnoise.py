import numpy as np
import argparse

#def gausian_ppm(arshape,ppm=500.0):
#    sigma=ppm*1.e-6
#    print(arshape)
#    noise=np.random.standard_normal(arshape)*sigma
#    return noise

def ppm2nphoton(ppm):
    return 1.0/(ppm*1.e-6)**2

def convfactor_photon(ppm,meanar):
    cfac=ppm2nphoton(ppm)/meanar
    return cfac

