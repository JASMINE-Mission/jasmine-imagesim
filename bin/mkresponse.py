#!/usr/bin/env python3.8
"""Make a response distribution

  usage:
    mkresponse.py [-h|--help] -t tel -s src -d det -r res

 options:
   --help      show this help message and exit
   -t tel      Telescope parameter
   -s src      Source (Object) parameter
   -d det      Detectir parameter
   -r res      output Responce file
"""
from docopt import docopt             # command line interface
import json
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
  k,WLdefined,EPdefined,WLshort,WLlong=extract_json.exttel(tel)
  WLdet,QEdet=extract_json.extQE(det)
  Rv, JH, alp=extract_json.extsrc(src)
  
  #calc response
  Tr,WL,Npr=response.calc_response(Rv, JH, alp, k,\
                      WLdefined,EPdefined,WLshort,WLlong,\
                      WLdet,QEdet)

  #fits
  with open(args['-r'],mode="w") as f:
    f.write('{\n')
    f.write('\"Ntot":{\n')
    f.write('  \"title\" : \"0 mag Total electrons/m2/sec \",\n')
    f.write('  \"val\" : {:.8e}'.format(Tr)+' },\n')
    f.write('\"WLdef\":{\n')
    f.write('  \"title\" : \"wavlength points\",\n')
    for i in range(len(WL)):
      f.write('  \"v{:02d}\" : {:.8e}'.format(i,WL[i]))
      if i < len(WL)-1 :
        f.write(',\n')
      else :
        f.write('},\n')
    f.write('\"SPR\":{\n')
    f.write('  \"title\" : \"Spectal distribution of input photons\",\n')
    for i in range(len(WL)):
      f.write('  \"v{:02d}\" : {:.8e}'.format(i,Npr[i]))
      if i < len(WL)-1 :
        f.write(',\n')
      else :
        f.write('}\n')
    f.write('}\n')
