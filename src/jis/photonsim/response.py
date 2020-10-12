import numpy as np
import sys
import math
from scipy import interpolate
from scipy import integrate   

# physical constants
c = 299792458.     ; # [m/s] 
h = 6.62607015e-34 ; # [J/s] 
# MKO NIR system from 2005PASP..117..421T, 2005PASP..117.1459T
# wavelength [um]
WL_MKO   = np.array([ 0.5450 , 1.250  , 1.644  , 2.121  , 2.149  , 2.198  ,
                      3.754  , 4.702])
# flux [W/m2/um]
FLUX_MKO = np.array([3.68e-08,3.01e-09,1.18e-09,4.57e-10,4.35e-10,4.00e-10,
                     5.31e-11,2.22e-11])
# flux [Photons/m2/um]
NP_MKO = FLUX_MKO * WL_MKO * 1.0e-6 / h / c ;
# interpolate flux
logL = np.log10(WL_MKO)
logN = np.log10(NP_MKO)
Npspline = interpolate.interp1d(logL,logN,kind="cubic")
# A(lambda)/Av
AWL = lambda wl , R : EXT_a(wl) + EXT_b(wl)/R
# Band definition
WL_J = 1.250
WL_H = 1.644

def calc_response(Rv, JH, alp, k, WLdefined,EPdefined,WLshort,WLlong,WLdet,QEdet):
    # from Rv and J-H, calculate Av
    JHA = AWL(WL_J,Rv) - AWL(WL_H,Rv)
    Av = JH/JHA
    
    # Array for wavelength
    # WL[i] should be in (WLshort,WLlong and 0.1 um step)
    Wlist=[]
    eps = 1e-8 
    for i in range(22):
        w = i/10+0.4  # from 0.4 um to 2.5 um
        if i/10+0.4 >= WLshort-eps and i/10+0.4 <= WLlong+eps :
            Wlist.append(i/10+0.4)
    WL = np.array(Wlist)

    # ingerpolate efficiency and qe
    EPinter = interpolate.interp1d(WLdefined,EPdefined,kind='linear')
    QEinter =  interpolate.interp1d(WLdet,QEdet,kind='linear')
    EP = EPinter(WL)
    QE = QEinter(WL)
    # 0 mag Photon Flux
    #print(WL)
    # Photon Fulx
    Np=Nphotons(WL)
    # band definition
    NpJ=Nphotons(WL_J)
    NpH=Nphotons(WL_H)
    NpHw = NpJ*alp+NpH*(1-alp) # Hw band  0 mag flux
    NprJ = NpJ*math.pow(10.0,-AWL(WL_J,Rv)*Av/2.5) # J band reddened flux
    NprH = NpH*math.pow(10.0,-AWL(WL_H,Rv)*Av/2.5) # H band reddened flux
    NprHw = NprJ*0.75+NprH*0.25 # Hw band reddened flux
    # reddenend number of photons
    Npr = np.empty(len(WL))
    for i in range(len(WL)):
        Npr[i] = EP[i]*QE[i]*Np[i]*math.pow(10.0,-AWL(WL[i],Rv)*Av/2.5)*NpHw/NprHw
    # Total Photons
    Tr = integrate.simps(Npr,WL)
    return Tr,WL,Npr

def Nphotons(WL):   # interpolated flux from MKO-NIR system
  logwl = np.log10(WL)
  lognp = Npspline(logwl)
  npw = np.power(10,lognp)
  return npw

# Interstellar extinction curve coefficients from 1989ApJ...345..245C
def EXT_a(wl):
  x = 1/wl
  if x >= 0.3 and x<=1.1 :
    return 0.574 * math.pow(x,1.61)
  elif x <= 8. :
    y = x - 1.82
    y1 = 1+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4
    y2 =   0.01979*y**5-0.77530*y**6+0.32999*y**7
    return y1 + y2 
  else :
    print('WL out of range error')
    sys.exit()
    
def EXT_b(wl):
  x = 1/wl
  if x >= 0.3 and x<=1.1:
    return -0.527* math.pow(x,1.61)
  elif x <= 8.:
    y = x - 1.82
    y1 = 1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5
    y2 = 5.30260*y**6-2.09002*y**7
    return y1 + y2
  else :
    print('WL out of range error')
    sys.exit()
