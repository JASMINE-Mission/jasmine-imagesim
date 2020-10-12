import numpy as np
import argparse
import pylab
import matplotlib.pyplot as plt

    
def gentraj_k2like(ntime,basepos,nsub,basesig=0.1,lpix=1.0,lpixsig=0.1,pixret=30,grad=0.74,gradsig=0.1,seed=None):
    nret=np.int(float(ntime)/float(pixret*nsub))+1
    np.random.seed(seed)
    karray=np.random.randn(nret)*gradsig+grad
    larray=np.random.randn(nret)*lpixsig+lpix
    bx=basepos[0]+np.random.randn(nret)*basesig
    by=basepos[1]+np.random.randn(nret)*basesig
    
    thetax=np.array([])
    thetay=np.array([])
    for i in range(0,nret):
        k=karray[i]
        l=larray[i]
        deltax=l/np.sqrt(1.0+k*k)
        deltay=k*deltax
        tx=np.linspace(bx[i],bx[i]+deltax,pixret*nsub).astype(np.float32)
        ty=np.linspace(by[i],by[i]+deltay,pixret*nsub).astype(np.float32)
        thetax=np.hstack([thetax,tx])
        thetay=np.hstack([thetay,ty])
    theta=np.array([thetax[0:ntime],thetay[0:ntime]])
    return theta

def gentraj_encircled_k2like(ntime, basepos, nsub, radius=1.0, lpix=1.0, lpixsig=0.1, pixret=30, \
                             grad=0.74, gradsig=0.01,seed=None):
    nret = np.int(float(ntime) / float(pixret * nsub)) + 1
    karray = np.random.randn(nret) * gradsig + grad
    larray = np.random.randn(nret) * lpixsig + lpix

    thetax = np.array([])
    thetay = np.array([])
    j=0
    np.random.seed(seed)
    print("radius=",radius)
    while j<=ntime:
        for i in range(0, nret):
            k = karray[i]
            l = larray[i]
            deltax = l / np.sqrt(1.0 + k * k)
            deltay = k * deltax
            rx, ry = 2.0 * (np.random.rand(2) - 0.5)*radius
            tx = np.linspace(rx, rx+deltax, pixret*nsub).astype(np.float32)
            ty = np.linspace(ry, ry+deltay, pixret*nsub).astype(np.float32)
            mask=(tx*tx+ty*ty<radius*radius)
            if len(tx[mask]) > 1:
                thetax = np.hstack([thetax, tx[mask] + basepos[0]])
                thetay = np.hstack([thetay, ty[mask] + basepos[1]])
                j = j + len(tx[mask])

    theta = np.array([thetax[0:ntime], thetay[0:ntime]])
    return theta


def gentraj_random(ntime,basepos,nsub,basesig=1.0,subsig=0.1,seed=None):
    ndata=np.int(ntime/nsub)
    np.random.seed(seed)
    thetax = np.array([])
    thetay = np.array([])
    retposx=np.random.randn(ndata)*basesig
    retposy=np.random.randn(ndata)*basesig
    for i in range(0,ndata):
        thetax=np.hstack([thetax,basepos[0]+retposx[i]+np.random.randn(nsub)*subsig])
        thetay=np.hstack([thetay,basepos[1]+retposy[i]+np.random.randn(nsub)*subsig])

    theta=np.array([thetax[0:ntime],thetay[0:ntime]])
    return theta


############################################################
def gentraj_random_old(ntime,basepos,basesig=1.0,seed=None):
    np.random.seed(seed)
    thetax=basepos[0]+np.random.randn(ntime)*basesig
    thetay=basepos[1]+np.random.randn(ntime)*basesig
    theta=np.array([thetax[0:ntime],thetay[0:ntime]])
    return theta

def gentraj_ranwalk_old(ntime,basepos,dr=0.1,seed=None):
    np.random.seed(seed)
    phi=np.random.rand(ntime)*2*np.pi
    dx=dr*np.cos(phi)
    dy=dr*np.sin(phi)
    theta=np.array([np.cumsum(dx)+basepos[0],np.cumsum(dy)+basepos[1]])
    return theta

def gentraj_encircled_random_old(ntime,basepos,radius=1.0,seed=None):
    thetax=[]
    thetay=[]
    j=0
    np.random.seed(seed)
    while j<ntime:
        rx,ry=2.0*(np.random.rand(2)-0.5)
        rr=rx*rx+ry*ry
        if rr<1.0:
            thetax.append(rx * radius + basepos[0])
            thetay.append(ry * radius + basepos[1])
            j=j+1

    theta = np.array([thetax[0:ntime], thetay[0:ntime]])

    return theta

def gentraj_encircled_k2like_old(ntime, basepos, radius=1.0, lpix=1.0, lpixsig=0.1, pixret=30, \
                             grad=0.74, gradsig=0.01,seed=None):
    nret = np.int(float(ntime) / float(pixret)) + 1
    karray = np.random.randn(nret) * gradsig + grad
    larray = np.random.randn(nret) * lpixsig + lpix

    thetax = np.array([])
    thetay = np.array([])
    j=0
    np.random.seed(seed)
    print("radius=",radius)
    while j<=ntime:
        for i in range(0, nret):
            k = karray[i]
            l = larray[i]
            deltax = l / np.sqrt(1.0 + k * k)
            deltay = k * deltax
            rx, ry = 2.0 * (np.random.rand(2) - 0.5)*radius
            tx = np.linspace(rx, rx+deltax, pixret).astype(np.float32)
            ty = np.linspace(ry, ry+deltay, pixret).astype(np.float32)
            mask=(tx*tx+ty*ty<radius*radius)
            if len(tx[mask]) > 1:
                thetax = np.hstack([thetax, tx[mask] + basepos[0]])
                thetay = np.hstack([thetay, ty[mask] + basepos[1]])
                j = j + len(tx[mask])

    theta = np.array([thetax[0:ntime], thetay[0:ntime]])
    return theta

if __name__ == "__main__":

    
    #    theta=gentraj_encircled_random(1024, [6,6], radius=1.0)
    nsub=16
    ntime=1024*nsub
    basepos=[6, 6]
    #theta=gentraj_k2like(ntime,basepos,nsub)
    #theta=gentraj_random(ntime,basepos,nsub)
    pixsec=29.0 #arcsec/pix
    sigsec_interframe=0.0
    sigsec_intraframe=11.0 #arcsec

    theta=gentraj_random(ntime,basepos,nsub,basesig=sigsec_interframe/pixsec,subsig=sigsec_intraframe/pixsec,seed=None)
#    theta = gentraj_encircled_k2like(ntime,basepos,nsub, radius=2.0, lpix=1.4,gradsig=3.0)
    fig=plt.figure(figsize=(10,10))
    plt.gca().invert_yaxis()
    plt.plot(theta[0,::], theta[1,:], ".", c="green", markersize=2)
    plt.plot(theta[0,::nsub],theta[1,::nsub],"+",c="red",lw=1)
    plt.show()

    
    #theta=gentraj_k2like(ntime,basepos,basesig=0.1,pixret=30,lpix=1.0,lpixsig=0.1,grad=0.74,gradsig=0.01)
