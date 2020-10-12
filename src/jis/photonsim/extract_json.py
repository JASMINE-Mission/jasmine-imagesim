import numpy as np

def extsp(sp):
    k=len(sp['WLdef'])-1
    WL = np.empty(k)
    for i in range(k):
        WL[i] = sp['WLdef']['v{:02d}'.format(i)]
    NP = np.empty(k)
    for i in range(k):
        NP[i] = sp['SPR']['v{:02d}'.format(i)]
    Ntot = sp['Ntot']['val']
    return k,WL,NP,Ntot

def exttel(tel):
    # from Eopt in tel, prepare arrays for WL and E
    k = int((len(tel['Eopt'])-1)/2 ) # number of definition points
    WLdefined = np.empty(k) # Wave length
    EPdefined = np.empty(k) # Efficiency of optics
    for i in range(k):
        WLdefined[i] = tel['Eopt']['W{:02d}'.format(i)]
        EPdefined[i] = tel['Eopt']['V{:02d}'.format(i)]
    WLshort = WLdefined[0]
    WLlong =  WLdefined[k-1]
    return k,WLdefined,EPdefined,WLshort,WLlong

def extQE(det):
    # read QE from det
    k = int((len(det['QE'])-1)/2 ) # number of definition points
    WLdet = np.empty(k) # Wave length
    QEdet = np.empty(k) # Efficiency of optics
    for i in range(k):
        WLdet[i] = det['QE']['W{:02d}'.format(i)]
        QEdet[i] = det['QE']['V{:02d}'.format(i)]
    return WLdet, QEdet

def extsrc(src):
    Rv = src['Rv']['val']
    JH = src['J-H']['val']
    alp = src['Hwband']['val']
    return Rv, JH, alp
