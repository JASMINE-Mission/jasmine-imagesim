import math
import numpy as np
import pyfftw
import sys
import matplotlib.pyplot as plt

def calc_ace(rg,N,T,ace):
  """
  Simulate one dimensional attitude control error data
  """
  dT = T/N # time step
  # Set parameters and PSD function
  title = ace['title']
  # print("used model {0},{1},{2}".format(title,ace['model'],ace['note']))
  if title == 'ace001':
    f0  = ace['f0']
    alp = ace['alp']
    bet = ace['bet']
    nmax = ace['n']
    P = np.empty(nmax)
    F = np.empty(nmax)
    H = np.empty(nmax)
    for i in range(nmax):
      pi = "p{}".format(i+1)
      hi = "h{}".format(i+1)
      fi = "f{}".format(i+1)
      P[i] = ace[pi]
      H[i] = ace[hi]
      F[i] = ace[fi]
  else :
    print("Not supported title {0}" .format(title))
    sys.exit()

  if title == 'ace001' :
    """
    Set Fourier domain ACE

    ACE(t) should be a real (not complex) function.
    Fourier domain ACE F(f) should...
      F(0 or N/2 ) should be real
      F(N-f) = conjugete( F(f) )
        then 
          ReF(-f) = ReF(f) ,  ImF(-f) = -ImF(f)
    """
    # initialize data for fft
    data = pyfftw.zeros_aligned((N),dtype='complex128')
    # independent data
    for t in range(int(N/2)+1):
      f = t/(N*dT)
      psd = 1/f0/(1+(f/f0)**2)
      for i in range(nmax):
        psd = psd + P[i]/H[i]/(1+((f-F[i])/H[i])**2)
      s0 =  math.sqrt(psd)
      sig = math.pow(s0,bet)
      noise = rg.normal(scale=sig)
      s = s0 + alp * noise
      if (t==0 or t==N/2) : # real
        th = 0.0
      else :
        th = rg.uniform() * 2*math.pi
      data.real[t] = s*math.cos(th)
      data.imag[t] = s*math.sin(th)

    # dependent data
    for t in range(int(N/2)-1):
      tt = N-t-1
      data.real[tt] =  data.real[t+1]
      data.imag[tt] = -data.imag[t+1]

    # iFFT
    ft = pyfftw.interfaces.numpy_fft.ifft(data)
    # statistics
    mean  = np.mean( ft.real )
    std = np.std( ft.real, ddof=1) # sqrt( SUM((x-<x>)**2)/(N-1) )
    acedata = ( ft.real - mean ) / std
    psdn = data.real*data.real+data.imag*data.imag

    return acedata,psdn

def plot_ace(N,T,data,psdn,plotfile):
  dT = T/N # time step
  fig = plt.figure()
  ax1 = fig.add_subplot(2,1,1)
  Ta = np.arange(0,N*dT,dT)
  p1 = ax1.plot(Ta,data)
  sig = np.std( data, ddof=1)
  ax1.set_ylabel('Atitude Control Error [arcsec]')
  ax1.set_xlabel('Time [sec]')
  ax1.text(dT*20,3.0,"Sigma={:.2f}".format(sig))
  ax1.set_ylim(bottom=-4,top=4)
  ax2=fig.add_subplot(2,1,2)
  FR= np.arange(0,int(N/2+1))
  FR= FR /(N*dT)
  p2=ax2.plot(FR,psdn[0:int(N/2+1)],marker='o',markersize=1,linestyle='None')
  ax2.set_xscale('log')
  ax2.set_yscale('log')
  ax2.set_ylabel('PSD')
  ax2.set_xlabel('Frequency [Hz]')

  plt.subplots_adjust(hspace=0.4)

  fig.savefig(plotfile,transparent=True)

  plt.show()


def ace2d(x_scale,y_scale,pix_scale,N,xdata,ydata,xN):
  """
  Make an atitude control error map
  
  xdatar : X-axis simulated ACE array
  ydata  : Y-axis simulated ACE array
  xN     : simulated ACE array size
  x_scale :  scale of xdata 
  y_scale :  scale of ydata 
  pix_scale : pixel scale of Output image
  N : output image is N x N 
  """
  data = np.zeros( (N,N) , dtype=int)
  for i in range(xN):
    x = xdata[i]*x_scale/pix_scale
    y = ydata[i]*y_scale/pix_scale
    xi = math.floor(x + 0.5 + N/2 )
    yi = math.floor(y + 0.5 + N/2 )
    data[yi,xi] = data[yi,xi]+1

  return data
