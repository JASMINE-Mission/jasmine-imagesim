#!/usr/bin/env python3.8
"""Make a response distribution

.. code-block:: bash

  usage:
    mkresponse.py [-h|--help] -t tel -s src -d det -r res

  options:
    --help  show this help message and exit
    -t tel  Telescope parameter
    -s src  Source (Object) parameter
    -d det  Detectir parameter
    -r res  output Responce file
"""
from docopt import docopt             # command line interface
import json
import numpy as np
from jis import photonsim
from jis.photonsim import response
from jis.photonsim import extract_json


#  command line interface
if __name__ == '__main__':
  args = docopt(__doc__)
  # get parameter
  with open(args['-t']) as f:
    tel = json.load(f)
  with open(args['-s']) as f:
    src = json.load(f)
  with open(args['-d']) as f:
    det = json.load(f)

  #extract values from json/array
  WLdefined,EPdefined=extract_json.extTel(tel)
  WLdet,QEdet=extract_json.extQE(det)
  Rv, JH, alp=extract_json.extSrc(src)

  #calc response
  Tr,WL,Npr=response.calc_response(
    Rv,JH,alp,WLdefined,EPdefined,np.min(WLdefined),np.max(WLdefined),WLdet,QEdet)

  data = {
    "Ntot": {
      "title": "Zero-mag total detected photon (electron) rate (e-/m^2/sec).",
      "val": Tr,
    },
    "wavelength": {
      "title": "Wavelength grid in micron.",
      "val": list(WL),
    },
    "spectrum": {
      "title": "Spectal distribution of detected photon (electron) rate (e-/m^2/sec/um).",
      "val": list(Npr),
    },
  }

  #fits
  with open(args['-r'],mode="w") as f:
    f.write(json.dumps(data, indent=True))
