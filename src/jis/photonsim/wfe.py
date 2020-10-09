import math
import numpy as np
from jis import photonsim
from jis.photonsim import zernike                       # Zernike polynomials functions
import json

def calc_wfe(EPD,efile):
    # make a little bit larger data
    N = int(EPD + 4) 
    data = np.zeros((N,N),dtype=float)
    
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
            
            
    for i in range(NP):
        for iy in range(N):
            for ix in range(N):
                y = iy - N/2
                x = ix - N/2
                rho = math.sqrt( y*y+x*x ) / (EPD/2) 
                if rho <= 1 :
                    th = math.atan2( y, x)
                    data[iy,ix] = data[iy,ix] + Za[i]*zernike.Zernike(Zn[i],Zm[i],rho,th)
                else :
                    data[iy,ix] = math.nan

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
