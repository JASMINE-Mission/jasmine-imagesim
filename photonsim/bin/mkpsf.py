#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""
  Make a PSF

  usage:
    mkpsf.py [-h|--help] -a apt -w wfe -s spec -c cntl -p psf.fits

 options:
   --help         show this help message and exit
   -a apt         aperture mask
   -w wfe         wfe map or null
   -s spec        electrons/sec/um at lambada
   -c cntl        control parameter
   -p psf.fits
"""
from docopt import docopt             # command line interface
import math
import numpy as np
import astropy.io.fits as fits
import json
import pyfftw
import sys

#  command line interface
if __name__ == '__main__':
  args = docopt(__doc__)

# spectral response
with open(args['-s']) as f:
  sp = json.load(f)

k=len(sp['WLdef'])-1
WL = np.empty(k)
for i in range(k):
  WL[i] = sp['WLdef']['v{:02d}'.format(i)]
NP = np.empty(k)
for i in range(k):
  NP[i] = sp['SPR']['v{:02d}'.format(i)]
Ntot = sp['Ntot']['val']


# get the aperture mask
ahdul = fits.open(args['-a'])
ahdr = ahdul[0].header
aN = ahdr['NAXIS1']
if ahdr['NAXIS1'] != ahdr['NAXIS2']:
  print("NAXIS1 {0} != NAXIS2 {1}".format(aN,ahdr['NAXIS2']))
  sys.exit()
Stel = ahdr['STEL']   # Area of collecting photons
adata = ahdul[0].data

# get the WFE map
if args['-w'] == 'null':
  wN = aN
  wfe = np.zeros((wN,wN))
else:
  whdul = fits.open(args['-w'])
  whdr = whdul[0].header
  wN = whdr['NAXIS1']
  if whdr['NAXIS'] != 2:
    print("NAXIS = {0} is bad in {1}".format(whdr['NAXIS'],args['-e']))
    sys.exit(-1)
  if whdr['NAXIS1'] != whdr['NAXIS2']:
    print("NAXIS1 {0} != NAXIS2 {1}".format(wN,whdr['NAXIS2']))
    sys.exit(-1)
  wfe = whdul[0].data *2*math.pi

wfer = np.nan_to_num( wfe )
wfec = np.empty( (wN,wN),dtype='complex128')
# 波長を k 個の点で規定している。
# 口径 D で、計算領域の大きさを N とするとき、
# フーリエ変換で得られる PSF は D/N x lambda/D=lambda/N の角度スケールになる。 
# 異なるいくつかの波長 WL0,WL1,,,WLn で生成したPSFのセルスケールを
# 合わせようとおもったら、計算領域を波長に比例させて
# N0 = M WL0 , N1 = M WL1 ,,, Nn = M WLn 
# と選べばよい。
# control parameter
with open(args['-c']) as f:
  cp = json.load(f)
M = cp['M']['val']
N0 = 520
image = np.zeros( (N0,N0) ,dtype ='float' )
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
  ft = pyfftw.interfaces.numpy_fft.fft2(data)
  fts = np.fft.fftshift(ft)

  # 結果を入射光子数(電子数)の重みを付けて加算する。
  i1 = int(N/2-N0/2)
  i2 = int(N/2+N0/2)
  NPm = (NP[i] + NP[i+1])/2
  image = image + NPm*(fts.real[i1:i2,i1:i2]**2 + fts.imag[i1:i2,i1:i2]**2)
 
# 最後に規格化しておく、値は1秒あたりのelectron数になるように
s = np.sum(image)
image = image/s * Ntot * Stel
# Save 
hdu = fits.PrimaryHDU(image)
hdu.header['NTOT'] = Ntot  # total flux/m2/sec

hdu.header['APTFILE'] = ahdr['APTFILE']     # telescope.json
hdu.header['EPD']     = ahdr['EPD']         # pupil diameter (mm)
hdu.header['COBS']    = ahdr['COBS']        # central obscuration ratio 
hdu.header['STYPE']   = ahdr['STYPE']       # Spider type
hdu.header['TSP']     = ahdr['TSP']         # spider thickness (mm)
hdu.header['STEL']    = ahdr['STEL']        # Area of telescope aperture (m2)
hdu.header['SPRES']   = args['-s']          # spectral responce 
hdu.header['M'] = M   # Number of FFT cells per wavelength in um

hdulist = fits.HDUList([hdu])
hdulist.writeto(args['-p'],overwrite=True)
