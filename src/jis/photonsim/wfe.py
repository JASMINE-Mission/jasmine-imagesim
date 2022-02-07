import os, glob
import math
import numpy as np
import jis
from jis import photonsim
from jis.photonsim import zernike                       # Zernike polynomials functions
import json
import pandas as pd
from scipy.interpolate import RectBivariateSpline

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


def calc_wfe_fringe37(EPD, filename, scale, positions, omit_tilt=True):
    """
    This function calculates wavefront error (wfe)
    at each position in the positions array. The wfe is
    calculated by taking summation of zernike polynomials.
    The amplitudes are given by the file whose name is
    given by the filename parameter. The file should have

    The file should be a csv file containing
    2D Zernike coefficient data. The first column
    of the csv data must be 'xan', 'yan', 1, 2, ..., 37.
    'xan' and 'yan' represent positions on the focal plane in deg.
    1, 2, ..., 37 are the Fringe Zernike indices.
    Other columns must have the positions (row1, row2) and
    Zernike coefficients (row3 to row39; to be scaled
    by the scale parameter).

    Args:
        EPD         (float): Entrance pupil diameter (apt-cell=mm).
        filename      (str): Filename of the csv file having the Zernike coefficients.
        scale       (float): Scaling factor in um. This will be multiplied to the Zernike coefficients.
        positions (ndarray): Array of positions ([[x0, y0], [x1, y1], ..., [xn, yn]]) in deg.
        omit_tilt    (bool): If true, tilt terms are omitted (default: True).
        *** Apt-cell scale is assumed to be 1 mm/apt-cell!!! ***
        *** Coefficients are assumed to be in um!!!  ***

    Returns:
        wfe (ndarray): (# of positions) x EPD+4 x EPD+4 data array of the calculated wfe.
                       The unit is um.

    """

    # Making a dictionary of functions which calculate
    # amplitudes at each position.
    dname = os.path.dirname(jis.__file__)
    fname = os.path.join(dname, filename)
    if glob.glob(fname):
        za_functions = read_FringeZernike37(fname, scale)    # Searching in jis dir.
    else:
        za_functions = read_FringeZernike37(filename, scale) # Searching in other dir.

    # Setting (j, n, m) indices in the Fringe37 convention.
    indices = zernike.FringeID37()

    if omit_tilt:
        indices = indices[indices['j']>=4]
    else:
        indices = indices[indices['j']>=2]

    # Calculating wfe.
    wfe = []
    for p in positions:
        Za = []
        for j in indices['j']:
            Za.append(za_functions[j](p[0], p[1])[0])
        Za = np.array(Za)
        wfe.append(calc_wfe_from_Zernike_param_array(EPD, indices['n'], indices['m'], Za))
    wfe = np.array(wfe)

    return wfe


def calc_wfe_from_Zernike_param_array(EPD, Zn, Zm, Za):
    """
    This function calculates wavefront error map from
    given Zernike parameters.

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


def read_FringeZernike37(filename, scale):
    '''
    This function reads a csv file containing
    2D Zernike coefficient data and returns a
    dictionary of interpolation functions.

    The first column of the csv data must be
    'xan', 'yan', 1, 2, ..., 37. 'xan' and 'yan'
    represent positions on the focal plane in deg.
    1, 2, ..., 37 are the Fringe Zernike indices.
    Other columns must have the positions (row1, row2) and
    Zernike coefficients (row3 to row39; to be scaled
    by the scale parameter).

    The returned dictionary has 37 functions.
    The keys are 1, 2, ..., 37 which correspond to the
    Fringe Zernike indices. The function that is stored
    as j is the interpolation function of the
    j-th Fringe Zernike coefficient.
    If you want to have the j-th coefficient
    at a position of (x, y) (in deg), it is calculated
    by functions[j](x, y).

    Args:
        filename (str): Filename of the input csv file.
        scale (float) : Scaling factor (um).

    Returns:
        functions (dict): Dictionary of the interpolation functions.

    '''

    # Loading data.
    df = pd.read_csv(filename, index_col=0, header=None).T

    # Load grid data.
    xan = np.array(df['xan'])
    yan = np.array(df['yan'])
    # RectBivariateSpline requires regular grid points.
    # The index `idx` is used to rearrange the Zernike coefficients
    # to make the functions compatible with the `interp2d` results.
    xtic = np.sort(np.unique(xan))
    ytic = np.sort(np.unique(yan))
    idx = np.argsort(1e3*xan+yan)
    ny,nx = ytic.size,xtic.size

    # Making and storing interpolation functions.
    functions = {}
    for i in range(1,38):
        zarr = np.array(df[str(i)])[idx].reshape((nx,ny)) * scale
        functions[i] = RectBivariateSpline(xtic, ytic, zarr, kx=1, ky=1)

    return functions
