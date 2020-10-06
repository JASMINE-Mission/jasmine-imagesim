import numpy as np
import sys

def dQ_const(x,rho,Ei0):
    #x = T/tau
    return Ei0*rho*(1.0 + (-1.0 + np.exp(-x))/x)

def persistence_const(x,rho,Ei0,Qij):
    Qij_prev=np.copy(Qij)    
    dQij=dQ_const(x,rho,Ei0)
    dQi=np.sum(dQij)
    Qij=dQij+Qij_prev*np.exp(-x)        
    Qi=np.sum(Qij)
    Ei = Ei0 - dQi + np.sum(Qij_prev*(1.0-np.exp(-x)))
    return Qij,Qi,Ei

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    #Reproduce Fig 23 of Tulloch+2019 1908.06469v1
    N=40
    T=600.0 #sec
    tframe=np.array(range(0,N))*T
    
    Qi=0.0
    E0=7.5e4*np.ones(N)
    tau=np.array([1.0,10.0,100.0,1000.0,10000.0])
    rho=np.array([1.e-3,1.5e-3,1.5e-3,2.e-3,3e-3])*0.8#*1.e2

    x=T/tau
    E=[]
    Q=[]
    Qij=np.zeros(len(tau)) #initial trapped charge

    #Up-the-ramp
    for i in range(0,N):
        Qij,Qi,Ei=persistence_const(x,rho,E0[i],Qij)
        E.append(Ei)
        Q.append(Qi)

    #Dark
    N=2000
    T=5.0 #sec
    x=T/tau

    tframe=np.concatenate([tframe,np.array(range(0,N))*T+tframe[-1]])
    Ed0=np.zeros(N)
    E0=np.concatenate([E0,Ed0])
    for i in range(0,N):
        Qij,Qi,Ei=persistence_const(x,rho,Ed0[i],Qij)
        E.append(Ei)
        Q.append(Qi)

    fs=18
    fig=plt.figure(figsize=(10,7))
    ax=fig.add_subplot(212)
    ax.plot(tframe,Q,".")
    plt.ylabel("trapped signal [e-]",fontsize=fs)
    plt.xlabel("time [sec]",fontsize=fs)
    plt.ylim(0,400.0)
    plt.xlim(0,35000.0)
    plt.tick_params(labelsize=fs)

    ax3=fig.add_subplot(211)
    ax3.plot(tframe,E0,label="input")
    ax3.plot(tframe,E,label="measured")
    plt.legend()
    plt.xlim(0,35000.0)
    plt.yscale("log")
    plt.ylabel("Input [e-]",fontsize=fs)
    plt.tick_params(labelsize=fs)
    #plt.title("Mocking Fig. 23 in Tulloch+2019")
    plt.savefig("tulloch_fig23_mock.png")
    plt.show()
