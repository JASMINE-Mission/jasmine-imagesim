#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""Simulate atitude control error

.. code-block:: bash

  usage:
    simace.py [-h|--help] -n N -t T -e ace.json -m ace.fits [-p plot.png]

  options:
    -h --help    show this help message and exit
    -n N         Length of output N steps
    -t T         total time
    -e ace.json  atitude control error parameter file
    -m ace.fits  output atitude control error
    -p plot.png  plot output

""" 
from docopt import docopt             # command line interface
import time
import numpy as np
import json
from jis.photonsim import ace
import astropy.io.fits as fits

# Random generator
#  ここでは、乱数の初期化は実行形式のプログラム中で行います。
#  シミュレーションをコントロールしている、なるべく外側で初期化する
#  ことによって、乱数を用いる複数の関数を使った時でも統一的に
#  乱数を選べるようにしようという考えによるものです。
t = round(time.time()*10)                    #  make a seed from time
rg = np.random.Generator(np.random.PCG64(t)) # set up the seed
"""
See https://numpy.org/doc/1.18/reference/random/index.html 
By default, Generator uses bits provided by PCG64 which has better 
statistical properties than the legacy MT19937 
"""

#  Command line interface
if __name__ == '__main__':
  args = docopt(__doc__)
  
  # Get parameters from command line
  N = int(args['-n'])   # Total number of time steps
  T = float(args['-t']) # Total time
  
#  Command line interface
if __name__ == '__main__':
  args = docopt(__doc__)
  
  # Atitude control error parameters from json file
  #  このパラメータは将来的に複雑化すると思われるので、ここでは
  #  json.load で得られる辞書型を渡すようにします。実質はファイル名を
  #  渡すのと変わりありませんが、ファイルアクセスは実行形式のほうで
  #  やっておこうという考えです。
  with open(args['-e']) as f:
    acep = json.load(f)
  
  # one dimentional ace
  data,psdn = ace.calc_ace(rg,N,T,acep)
  
  # Save ACE map
  hdu = fits.PrimaryHDU(data)
  hdu.header["ACE-FILE"] = args['-e']
  hdu.header["ACE-TOTT"] = T
  hdulist = fits.HDUList([hdu])
  hdulist.writeto(args['-m'],overwrite=True)
  
  # plot
  if args['-p'] :
    ace.plot_ace(N,T,data,psdn,args['-p'])
  
