"""
 Zernike polynomial code
"""
import math
import numpy as np
import sys

def Radial(n,m,rho):
    """
    This function returns Zernike amplitude at rho with index n and m.

    Zernike function is:
      Z^m_n(rho, phi) = R^m_n(rho) x cos(m phi)   (m>=0)
                        R^m_n(rho) x sin(|m| phi) (m<0),

    This function returns R^m_n(rho)

    Args:
        n, m (int)  : Indices of the Zernike term (n>=0; m<=n).
        rho  (float/ndarray): Distance from the center.

    Returns:
        r (float/ndarray): R^m_n(rho) of the Zernike term.

    Example:
        import numpy as np
        from jis.photonsim.zernike import radial

        zernike_r_2_0 = radial(2, 0, 0.3) 

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
            r = r + s*a/math.factorial(k)/b/c*rho**(n-2*k)
            s = -s 

    return r

def Zernike(n,m,rho,theta):
    """
    This function returns Zernike function values at (rho, theta)
    with indices of n and m, Z^m_n(rho, theta).

    Args:
        n, m  (int)  : Indices of the Zernike function.
        rho   (float/ndarray): Distance from the center.
        theta (float/ndarray): Azimuthal angle (rad).   

    Returns:
        z (float/ndarray): Zernike (m, n) value at (rho, theta).

    Example:
        from jis.photonsim.zernike import Zernike 

        z_2_0 = Zernike(2, 0, 0.3, 0.1)

    """
    mm = abs(m)
    if n < 0 or mm > n :
        print("Error n={} m={}".format(n,m))
        sys.exit()

    if m>=0 :
        z = Radial(n,mm,rho) * np.cos( m*theta )
    else :
        z = Radial(n,mm,rho) * np.sin( mm*theta )

    return z

def ZernikeID(N):
    """
    This function returns first N Zernike-indices n and m as arrays.

    Args:
        N (int): Number of the indices you want to obtain.

    Return:
        ZIDn (ndarray): Array of n-indices.
        ZIDm (ndarray): Array of m-indices.

    Example:
        from jis.photonsim.zernike import ZernikeID

        ns, ms = ZernikeID(6) 
         
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

    return ZIDn, ZIDm


def FringeID37():
    '''
    This function returns the 37 Fringe's indices.
    The Zernike polynomials are defined with radial indices n
    and azimuthal indices m. The Fringe convention defines
    an indexing rule where an index j is defined with n and m as
        j = (1+(n+|m|)/2)^2 - 2|m| + |sgn(m)|*(1-sgn(m))/2.
    In some softwares, wavefront patterns are described with
    the first 37 coefficients in the Fringe convention (j=1-37).
    This function returns the 37 combinations of (j, n, m).

    Args:
        None

    Returns:
        indices (ndarray): The 37 Fringe's indices.
                           dtype=[('j', '<i4'), ('n', '<i4'), ('m', '<i4')]

    '''

    arr = []
    for n in range(0, 13):
        for m in range(-n, n+1, 2):
            j = (1+(n+np.abs(m))/2)**2-2*np.abs(m)+np.abs(np.sign(m))*(1-np.sign(m))/2 # Fringe indexing.
            arr.append([j, n, m])
    arr = np.array(list(map(tuple, arr)), dtype=[('j','<i4'),('n','<i4'),('m','<i4')])
    indices = np.sort(arr, order='j')[0:37] # Last term should be (37, 12, 0) ???

    return indices
