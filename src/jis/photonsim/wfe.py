import math
import numpy as np
from jis import photonsim
from jis.photonsim import zernike                       # Zernike polynomials functions
import json

def calc_wfe(EPD,efile):
    """
    This function calculates wavefront error pattern
    by calculating summation of zernike polynomials.
    Coefficients and orders of the zernike terms are
    given by the input json file, efile.
    In the simulation, the unit of the coefficients and
    the calculated wfe is assumed to be um.

    Args:
        EPD   (float): Entrance pupil diameter (apt-cell=mm).
        efile (str)  : Filename of the json file having wavefront error parameters.
        *** Apt-cell scale is assumed to be 1 mm/apt-cell!!! ***
        *** Coefficients are assumed to be in um!!!  ***

    Returns:
        wfe (ndarray): EPD+4 x EPD+4 data array of the calculated wavefront error pattern.
                       The unit is assumed to be um in the simulation.

    Example:
        from jis.photonsim.wfe import calc_wfe

        wfe_pattern = calc_wfe(500, "wfe.json").

    """

    # Get wavefront error parameters from json file
    with open(efile) as f:
        p = json.load(f)
        NP = int(p['N-polys'])
        Zn = np.empty(NP)
        Zm = np.empty(NP)
        Za = np.empty(NP)
        for i in range(NP):
            Zn[i] = int(p['z{:03d}-n'.format(i+1)])
            Zm[i] = int(p['z{:03d}-m'.format(i+1)])
            Za[i] = float(p['z{:03d}-a'.format(i+1)])

    wfe = calc_wfe_from_Zernike_param_array(EPD, Zn, Zm, Za)

    return wfe


def calc_wfe_from_Zernike_param_array(EPD, Zn, Zm, Za):
    """
    This function calculates wavefront error map from
    given Zernike parametes.

    Args:
        EPD  (float): Entrance pupil diameter (apt-cell=mm).
        Zn (ndarray): 1D array containing radial indices n of the Zernike polynomials (Z^m_n).
        Zm (ndarray): 1D array containing azimuthal indices m of the Zernike polynomials (Z^m_n).
        Za (ndarray): 1D array containing amplitudes of the Zernike polynomials in um.

    Returns:
        wfe (ndarray): (EPD+4) x (EPD+4) array of the calculated wavefront error map in um.

    """

    # Making a little bit large array.
    N   = int(EPD+4.)
    wfe = np.zeros((N, N), dtype=float)

    # Calculating x, y, rho, and theta.
    y, x  = np.indices((N, N)) - N/2
    rho   = np.sqrt(y*y+x*x)/(EPD/2)
    theta = np.arctan2(y, x)

    # Calculating wfe map.
    wfe[rho>1] = np.nan
    pos = np.where(rho<=1)
    for i in range(np.size(Za)):
        wfe[pos] += Za[i]*zernike.Zernike(Zn[i], Zm[i], rho[pos], theta[pos])

    return wfe


def calc_dummy_wfe(EPD,efile):
    """
    This function calculates a dummy wavefront error pattern.
    The arguments are the same as `calc_wfe` but not used except for `EPD`.
    This returns an array, where the inside of the entrance pupil is filled with zeros.

    Args:
        EPD   (float): Entrance pupil diameter (apt-cell=mm).
        efile (str)  : Filename of the json file having wavefront error parameters.
        *** Apt-cell scale is assumed to be 1 mm/apt-cell!!! ***
        *** Coefficients are assumed to be in um!!!  ***

    Returns:
        data (ndarray): EPD+4 x EPD+4 data array of the calculated wavefront error pattern.
                        The unit is assumed to be um in the simulation.
    """
    # make a little bit larger data
    N    = int(EPD + 4)
    data = np.zeros((N,N),dtype=float)

    iy,ix = np.indices((N, N))
    y,x = iy - N/2, ix - N/2
    rho   = np.sqrt(y*y+x*x)/(EPD/2)

    data[rho>1] = np.nan
    return data



def wfe_model_z(rg,nmax,wlen,zodd,zeven):
  """
  wave front error のために zernike で n=2..(nmax-1) の 強度データを作成する。
  長さのスケールに用いる波長を wlen とする。
  n が奇数のときは、強度を wlen/zodd とし、偶数の時は wlen/zeven とする。
  そして、 m=0 ではこの強度をそのまま用い、 m がゼロでなければ、
  この強度を正の m の項に に cos(θ) 倍、 m が負の項に sin(θ) 倍したものを
  用いることにする。ここで θ は 0 から 2π の一様乱数である。
  """
  n=int(nmax*(nmax+1)/2)
  ZIDn,ZIDm = zernike.ZernikeID(n)
  k=1
  wfe={}
  wfe['title'] = 'Zernike polynomials'
  wfe['N-polys'] = n-3
  for i in range(3,n):
    if ZIDn[i]%2 == 1 :
      a=wlen/zodd
    else :
     a=wlen/zeven
    if ZIDm[i]<0:
      continue
    elif ZIDm[i]==0:
      wfe['z{:03d}-n'.format(k)] = int(ZIDn[i] )
      wfe['z{:03d}-m'.format(k)] = int(ZIDm[i] )
      wfe['z{:03d}-a'.format(k)] = a
      k=k+1
    else :
      th = rg.uniform() * 2*math.pi
      wfe['z{:03d}-n'.format(k)] = int(ZIDn[i])
      wfe['z{:03d}-m'.format(k)] = int(ZIDm[i])
      wfe['z{:03d}-a'.format(k)] = a*math.cos(th)
      k=k+1
      wfe['z{:03d}-n'.format(k)] = int(ZIDn[i])
      wfe['z{:03d}-m'.format(k)] = int(-ZIDm[i])
      wfe['z{:03d}-a'.format(k)] = a*math.sin(th)
      k=k+1
  return wfe
