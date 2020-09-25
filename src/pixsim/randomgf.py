#!/usr/bin/python
import sys
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pylab 

def simple_gp(nsamp,alpha):
    kmin=10
    kmax=int(nsamp/10.0)
    k=np.array(range(nsamp),dtype=np.float)
    f=(k**-alpha)/float(nsamp)
    f[k<kmin]=0.0
    f[k>kmax]=0.0
    print(f)
    sig=np.sqrt(f/2)
    a=np.random.normal(loc=0.0, scale=1.0, size=nsamp)*sig
    b=np.random.normal(loc=0.0, scale=1.0, size=nsamp)*sig
    c=np.hstack((a+b*1j,0.0+0.0j))
    delta=np.fft.irfft(c)
    #cumd=np.cumsum(delta)
    cumd=delta*float(nsamp)
    #amplitude
    x=np.array(cumd)*4.0

    return x,f



if __name__ == "__main__":


    N=1000
    alpha=3.0
    x,f=simple_gp(N,alpha)
    y,f=simple_gp(N,alpha)
    
    fig=plt.figure()
    ax=fig.add_subplot(121)
    ax.plot(f)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("frequency (1/cyc)")
    plt.ylabel("power")
    
    ax=fig.add_subplot(122)
#    plt.plot(x,color="C0")
    plt.plot(x,y,".",color="C0")
    plt.plot(x,y,alpha=0.3,color="C0")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()
