import numpy as np

def extsp(sp):
    """
    Summary:
        This function extracts spectral information 
        from a decoded spectral data, 'sp'.

    Args:
        sp (dict): Decoded json data of the spectral information.

    Returns:
        k    (int)    : Number of output data points.
        WL   (ndarray): Wavelength data.
        NP   (ndarray): Spectral response data (Num. of photon?).
        Ntot (float?) : Total... what? 

    Example:
        import json
        from jis.photonsim.extract_json import extsp

        fp = open("***.json")
        sp = json.load(fp) 
        k, WL, NP, Ntot = extsp(sp)

    """

    k  = len(sp['WLdef'])-1

    WL = np.empty(k)
    for i in range(k):
        WL[i] = sp['WLdef']['v{:02d}'.format(i)]

    NP = np.empty(k)
    for i in range(k):
        NP[i] = sp['SPR']['v{:02d}'.format(i)]

    Ntot = sp['Ntot']['val']

    return k, WL, NP, Ntot


def exttel(tel):
    """
    Summary:
        This function extracts telescope information
        from a decoded telescope data, 'tel'.

    Args:
        tel (dict): Decoded json data of the telescope information.

    Returns:
        k         (int)    : Number of output data points.
        WLdefined (ndarray): Wavelength data.
        EPdefined (ndarray): Optical efficiency.
        WLshort   (float)  : Shortest wavelength.
        WLlong    (float)  : Longest wavelength.        

    Example:
        import json
        from jis.photonsim.extract_json import exttel

        fp  = open("***.json")
        tel = json.load(fp) 
        k, WL, EP, WLmin, WLmax = exttel(tel)

    """

    # from Eopt in tel, prepare arrays for WL and E
    k = int((len(tel['Eopt'])-1)/2 ) # number of definition points

    WLdefined = np.empty(k) # Wavelength
    EPdefined = np.empty(k) # Efficiency of optics
    for i in range(k):
        WLdefined[i] = tel['Eopt']['W{:02d}'.format(i)]
        EPdefined[i] = tel['Eopt']['V{:02d}'.format(i)]

    WLshort = WLdefined[0]
    WLlong =  WLdefined[k-1]

    return k, WLdefined, EPdefined, WLshort, WLlong


def extQE(det):
    """
    Summary:
        This function extracts quantum efficiency information
        from a decoded detector data, 'det'.

    Args:
        det (dict): Decoded json data of the detector information.

    Returns:
        WLdet (ndarray): Wavelength data.
        QEdet (ndarray): Quantum efficiency of the detector.

    Example:
        import json
        from jis.photonsim.extract_json import extQE

        fp  = open("***.json")
        det = json.load(fp) 
        WL, QE = extQE(det)

    """ 

    # read QE from det
    k = int((len(det['QE'])-1)/2 ) # number of definition points

    WLdet = np.empty(k) # Wave length
    QEdet = np.empty(k) # Efficiency of optics
    for i in range(k):
        WLdet[i] = det['QE']['W{:02d}'.format(i)]
        QEdet[i] = det['QE']['V{:02d}'.format(i)]

    return WLdet, QEdet


def extsrc(src):
    """
    Summary:
        This function extracts object source information
        from a decoded object source data, 'src'.

    Args:
        src (dict): Decoded json data of the object source information.

    Returns:
        Rv  (float): Extinction factor Rv(=Av/E(B-V)).
        JH  (float): J-H color of... what?
        alp (float): Hw-band magnitude of the source?.

    Example:
        import json
        from jis.photonsim.extract_json import extsrc

        fp  = open("***.json")
        src = json(fp) 
        Rv, JH, alp = extsrc(src) 

    """
    
    Rv  = src['Rv']['val']
    JH  = src['J-H']['val']
    alp = src['Hwband']['val']

    return Rv, JH, alp
