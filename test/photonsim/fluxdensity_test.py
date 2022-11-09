import pytest
import numpy as np
from jis.photonsim.fluxdensity import flux_rJHKs_byTeff, Hw_absmag

def sort_check(f_array, mag_array):
    
    boolsort_fd = (f_array[0] < f_array[1]) & \
            (f_array[1] < f_array[2])

    boolsort_mg = (mag_array[0] > mag_array[1]) &\
            (mag_array[1] > mag_array[2])

    boolsort    = np.hstack((boolsort_mg, boolsort_fd))
    return np.all(boolsort)

def order_check(f_array):

    fmean       = np.mean(f_array, axis=1)
    flog_ar     = np.log10(fmean)

    b1  = (-13 < flog_ar[0]) & (flog_ar[0] < -12)
    b2  = (-10 < flog_ar[1]) & (flog_ar[1] < -9)
    b3  = (-9 < flog_ar[2]) & (flog_ar[2] < -8)

    return (b1 & b2) & (b2 & b3)

def test_fluxdensity():
    teff_list   = [3000.,6000.,9000.]

    f_array     = []
    mag_array   = []
    for t_eff in teff_list:
        w,f = flux_rJHKs_byTeff(t_eff)
        mag = Hw_absmag(t_eff)

        f_array.append(f)
        mag_array.append(mag)
    
    f_array     = np.array(f_array)
    mag_array   = np.array(mag_array)

    b1  = sort_check(f_array, mag_array)
    b2  = order_check(f_array)
    
    assert b1 & b2

if __name__=="__main__":
    print(test_fluxdensity())

    
