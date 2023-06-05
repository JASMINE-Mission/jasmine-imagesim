#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""Wave front error  random generattor

.. code-block:: bash

  usage:
    rnd_wfe.py [-h|--help] -n nmax -e z_even -o z_odd -z wfe.json

  options:
    -h --help    show this help message and exit
    -n nmax      Zernike max n
    -e z_even    Strength of Zernike poly. at even odder is 1.4/z_even
    -o z_odd     Strength of Zernike poly. at odd  odder is 1.4/z_odd 
    -z wfe.json  Generated Zernike polynomials strength

""" 
from docopt import docopt             # command line interface
import time
import numpy as np
import json
from jis.photonsim import zernike
from jis.photonsim import wfe

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
  nmax = int(args['-n'])   # Max order of Zernike polynomial
  zeve = float(args['-e']) 
  zodd = float(args['-o']) 
  wlen = 1.4
  
  wfe = wfe.wfe_model_z(rg,nmax,wlen,zodd,zeve)
  
  # Save 
  with open(args['-z'],mode='w') as f:
    json.dump(wfe, f, indent=2)
