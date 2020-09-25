import numpy as np
import argparse
import readflat as rf
import os
import simpix_stable as sp
import makeflat as mf
import matplotlib.pyplot as plt
import pylab
import transitmodel
import persistence
import sys
#
# X,Y in input/output are reverted in y, x in kernel
#
def dQ_const_array(x,rho,Ei0):
    #x = T/tau
    dim=np.shape(Ei0)
    return np.outer(Ei0.flatten(),rho*(1.0 + (-1.0 + np.exp(-x))/x)).reshape(dim[0],dim[1],len(x))


def persistence_const_array(x,rho,Ei0,Qij):
    Qij_prev=np.copy(Qij)    
    dQij=dQ_const_array(x,rho,Ei0)
    dQi=np.sum(dQij,axis=2)
    Qij=dQij+Qij_prev*np.exp(-x)

    Qi=np.sum(Qij,axis=2)
    Ei = Ei0 - dQi + np.sum(Qij_prev*(1.0-np.exp(-x)),axis=2)
    return Qij,Qi,Ei

def read_trapped_charge(Qtrap,ix,iy,pixdim,figsw=0,figtag=0):
    #get trapped charge in local view from ix,iy,pixdim
    sx=ix+np.int(pixdim[0])
    sy=iy+np.int(pixdim[1])
    if(sx<np.shape(Qtrap)[0] and sy<np.shape(Qtrap)[1]):
        trap_local=Qtrap[ix:sx,iy:sy,:]
    else:
        print("Invalid ix,iy")
        return None
    
    #######
    if figsw==1:
        a=plt.imshow(trap_local,interpolation="nearest",cmap="bwr")
        plt.colorbar(a)
        plt.title("sigma="+str(np.std(sflat)))
        plt.savefig("Qtrap"+str(figtag)+".png")
    #######

    return trap_local

def push_trapped_charge(Qtrap,ix,iy,pixdim,Qij):
    sx=ix+np.int(pixdim[0])
    sy=iy+np.int(pixdim[1])
    if(sx<np.shape(Qtrap)[0] and sy<np.shape(Qtrap)[1]):
        Qtrap[ix:sx,iy:sy,:]=Qij
    else:
        print("Invalid ix,iy")
        return None
    return Qtrap

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='This code generates the light curve with a planet assuming a random tracking in a detector plane with persistence. No interpix fluctuation. Intrapix model from CCD.')
    parser.add_argument('-t', nargs=1, default=["test.csv"],help='trajectory file', type=str)
    parser.add_argument('-s', nargs=1, default=[0.5],help='sigma of PSF gaussian. typical Kepler = 0.5 ', type=float)
    parser.add_argument('-f', help='flat correction?', action='store_true')
    args = parser.parse_args()
    sigma= args.s[0]
    sigma2=sigma*sigma
    trajfile=args.t[0]
    
    pixdim=[30,30] # dimension in the aperture pixel positions
    spixdim=[32,32] # subpixel dimension in a pixel
    ntime=1024

    #generate Gaussian flat 
#    flat = mf.gaussian_flat(sigma=0.01)
    flat = np.ones((1024,1024)) #UNIFORM FLAT
    
    gpixdim=np.shape(flat) # dimension for global pixel positions
    
    #intrapixel
    filex="intravx.csv"
    filey="intravy.csv"
    dirh="../data/intrapix"
    intrapix=rf.read_intrapix(filex,filey,spixdim,dirh)

    #trajectory
    sigmaXY=0.567
    theta_cen = np.array(pixdim)/2.0

    ndisturb=256
    nframe = 1500
    tframe = 7.0 # [sec] for 1 frame
    t=np.array(range(0,nframe))*tframe #sec

    #persistence model
    tau=np.array([1.0,10.0,100.0,1000.0,10000.0]) #sec
    rho=np.array([1.e-3,1.5e-3,1.5e-3,2.e-3,3e-3])*0.8*1.e2
    NQc = len(tau) #number of the defected pixel time scale
    Qtrap = np.zeros((gpixdim[0],gpixdim[1],NQc)) #trapped charge
    xtau=tframe/tau

    #planet model
    Rpin=1.0
    injlc, b=transitmodel.gentransit(t/3600/24,t0=75.0/60/24,Rp=Rpin)

    #initial global position
    x0=(0.5*(np.shape(flat)[0]-spixdim[0]))
    y0=(0.5*(np.shape(flat)[1]-spixdim[1]))
    x=x0
    y=y0
    
    lc=[]    
    theta = np.array([theta_cen])
    jx,jy=np.int(x),np.int(y)
    Qtsave=[]
    for iframe in range(0,nframe):

        #global position update
        dx,dy=np.random.rand(2)*0.01
#        dx,dy=0.0,0.0
        x,y=x+dx,y+dy

        jxp,jyp=jx,jy
        jx,jy=np.int(x),np.int(y)
        djx,djy=jx-jxp,jy-jyp
        interpix=rf.flat_interpix(flat,jx,jy,pixdim,figsw=0)

        #debug
        #jx,jy=np.random.rand(2)
        #jx=np.int(jx*(np.shape(flat)[0]-spixdim[0]))
        #jy=np.int(jy*(np.shape(flat)[1]-spixdim[1]))
        #interpix=rf.flat_interpix(flat,jx,jy,pixdim,figsw=0)
        
        #local position update
        theta = theta+np.array([[dx-djx,dy-djy]])
        theta_distrubed=  np.random.normal(0.0,sigmaXY,2*ndisturb).reshape((ndisturb,2)) + theta
        thetaT = theta_distrubed.T

        #debug
        #theta=  np.random.normal(0.0,sigmaXY,2*ndisturb).reshape((ndisturb,2)) + np.array([theta_cen])
        #thetaT = theta.T

        pixar=sp.simpix(thetaT,interpix,intrapix,sigma2=sigma2)

#        testimg=np.nansum(pixar,axis=2)
#        plt.figure()
#        plt.imshow(testimg)
#        plt.show()
#        sys.exit()
        
        #persistence
        Qij=read_trapped_charge(Qtrap,jx,jy,pixdim,figsw=0)        
        E0i=np.nansum(pixar,axis=2)
        Qij,Qi,Ei=persistence_const_array(xtau,rho,E0i,Qij)
        Qtrap = push_trapped_charge(Qtrap,jx,jy,pixdim,Qij)
        Qtsave.append(np.sum(Qtrap,axis=(0,1)))        
        lctmp=np.mean(np.sum(Ei))
        lc.append(lctmp)

    lc=lc/np.mean(lc)*injlc
    print("moving ",np.sqrt((x-x0)**2+(y-y0)**2),"pix")
    
    Qtsave=np.array(Qtsave)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    for ii in range(0,5):
        ax.plot(tframe*np.array(range(0,len(Qtsave[:,ii]))),Qtsave[:,ii],label="tau="+str(tau[ii])+" sec")
    plt.xscale("log")
    plt.xlabel("time [sec]")
    plt.ylabel("trapped charge")        
    plt.legend()
    plt.savefig("trapped_charge.png")
    plt.close()
    
    Nbin=43
    lcbin=[]
    tbin=[]
    for i in range(0,int(len(lc)/Nbin)):
        lcbin.append(np.mean(lc[i*Nbin:(i+1)*Nbin]))
        tbin.append(np.mean(t[i*Nbin:(i+1)*Nbin]))
    lcbin=np.array(lcbin)
    tbin=np.array(tbin)
        
    print("STD=",np.std(lcbin))

    ###
    np.savez("lc",[t/60,lc],[tbin/60,lcbin])
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(t/60,lc,".",alpha=0.3,label="frame")
    ax.plot(tbin/60,lcbin,"s",alpha=1.0,label="5 min binning")
    plt.legend()
    plt.xlabel("[min]")
    if args.f:
        plt.savefig("lc.flat.png")
    else:
        plt.savefig("lc.png")
