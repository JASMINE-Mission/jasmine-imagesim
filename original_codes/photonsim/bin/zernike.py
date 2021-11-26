"""
 Zernike polynomial code
"""
import math
import numpy as np
import sys

def Radial(n,m,rho):
  """
   Return Zernike amplitude at rho with index n and m
  """
  if n < 0 or m < 0 or m > n:
    print("Error n={} m={}".format(n,m))
    sys.exit()
  if (n-m)%2 != 0 :
    r = 0
  else :
    r = 0
    s = 1
    for k in range( int((n-m)/2) + 1 ) :
      a = math.factorial(int(n-k))
      b = math.factorial(int( (n+m)/2-k ))
      c = math.factorial(int( (n-m)/2-k ))
      r = r + s*a/math.factorial(k)/b/c*math.pow(rho,n-2*k)
      s = -s 
  return r

def Zernike(n,m,rho,theta):
  """
   Return Zernike value at (rho,theta) with index n and m
  """
  mm = abs(m)
  if n < 0 or mm > n :
    print("Error n={} m={}".format(n,m))
    sys.exit()
  if m>=0 :
    z = Radial(n,mm,rho) * math.cos( m*theta )
  else :
    z = Radial(n,mm,rho) * math.sin( mm*theta )
  return z

def ZernikeID(N):
  """
    Return first N Zernike-index n and m as arrays
  """
  ZIDn = np.empty(N+1,dtype=int)
  ZIDm = np.empty(N+1,dtype=int)
  id=0
  n = 0 
  m = 0
  ZIDn[id] = n
  ZIDm[id] = m
  while id<N  :
    if n == m :
      n=n+1
      m=-n
      id = id + 1
      ZIDn[id] = n
      ZIDm[id] = m
    else:
      m = m + 2
      id = id+1
      ZIDn[id] = n
      ZIDm[id] = m
  return ZIDn,ZIDm
