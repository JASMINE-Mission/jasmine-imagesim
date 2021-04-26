import numpy as np
import sys

def dQ_const(x,rho,Ei0):
    """
    Summary:
        This function calculates the increase of the trapped charge
        during an exposure time x (=Texp/tau), based on the model in Tulloch+19.
        This calculation assumes that pixel charge is zero at the beginning of the exposure.
        (See also, persistence_const and Kawahara-san's document.)

    Args:
        x   (ndarray): Exposure time normalized by the detrapping time constants.
                       The detector has various detrapping time constants (tau1, tau2,...).
                       If the exposure time is expressed as Texp,
                       x is described as [Texp/tau1, Texp/tau2,...].
        rho (ndarray): Fraction of trapped photocharge in each time constant bin at equilibrium.
        Ei0 (float)  : Photocharge on the pixel after the exposure.

    Return:
        dQ (ndarray): Increase of the trapped charge in each time constant bin.

    """
    dQ = Ei0*rho*(1.0 + (-1.0 + np.exp(-x))/x)

    return dQ

def persistence_const(x, rho, Ei0, Qij):
    """
    Summary:
        This function returns the trapped charges (Qij and Qi)
        and photocharge (Ei) after a time x (=t/tau).
        This calculation is based on the model by Tulloch+19.
        (for H2RG; See Tulloch+19 arXiv:1908.06469v1 and Kawahara-san's document.)

        For calculating the trapping process, please ensure that the photocharge
        is zero at the beggining of the exposure and set x and Ei0 to Texp/tau and
        the photocharge expected at the end of the exposure without trapping, respectively.
        Do not use this function to calculate an increase of the trapped charge
        in a small part of time in an exposure, because except for the first part,
        the photocharge is not zero (This is required to use dQ_const in the calculation).

        For calculating the detrapping process after a reset, you can divide the time
        after the reset and calculate the decrease of the trapped charge in each time bin
        by setting x and Ei0 to dt/tau and zero, respectively.

    Args:
        x   (ndarray) : A time normalized by the detrapping time constants.
                        The dector has various detrapping time constants (tau1, tau2...).
                        At a time t, x is described as [t/tau1, t/tau2,...].
                        For calculating the trapping process in an exposure time Texp,
                        set x to Texp/tau.
        rho (ndarray) : Fraction of trapped photocharge in each time constant bin at equilibrium.
        Ei0 (float)   : Photocharge on the pixel.
                        For calculating trapping process in an exposure,
                        set Ei0 to that expected in the end of the exposure without trapping.
                        For calculating detrapping process, set Ei0 to zero.
                        (If the photocharge is non-zero after a reset due to trapping,
                         does this function correctly work (TK)?)
        Qij (ndarray) : Initial trapped charge in each time constant bin.
                        "i" is the index of time. "j" is the index of the time constants.

    Returns:
        Qij (ndarray): Updated trapped charge in each time constant bin.
                       "i" is the index of time. "j" is the index of the time constants.
        Qi  (float)  : Updated total trapped charge.
        Ei  (float)  : Updated photocharge on the pixel.

    """

    Qij_prev = np.copy(Qij)           # Copying the initial trapped charge as Qij_prev.
    dQij = dQ_const(x, rho, Ei0)      # Calculating the increase of the trapped charge.
    dQi  = np.sum(dQij)               # Calculating the increase of the total trapped charge.
                                      # (\sum_{j} dQi_j)
    Qij  = dQij + Qij_prev*np.exp(-x) # Updated trapped charge in each time constant bin.
    Qi   = np.sum(Qij)                # Updated total trapped charge.
    Ei   = Ei0 - dQi + np.sum(Qij_prev*(1.0-np.exp(-x))) # Updated photocharge on the pixel.

    return Qij, Qi, Ei

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # Reproduce Fig 23 of Tulloch+2019 1908.06469v1
    N = 40    # Number of the up-the-ramp sequences.
    T = 600.0 # Exposure time of the sequence (sec).
    tframe = np.array(range(0,N)) * T

    Qi  = 0.0              # Initial trapped charge
    E0  = 7.5e4*np.ones(N) # Signal at the T-sec exp (75ke-/Tsec).
    tau = np.array([1.0,10.0,100.0,1000.0,10000.0])      # Detrapping time constants.
    rho = np.array([1.e-3,1.5e-3,1.5e-3,2.e-3,3e-3])*0.8
    # Fraction of the photocharge that is trapped in each time constant bin.

    x = T/tau
    E = []
    Q = []
    Qij = np.zeros(len(tau)) # Initial trapped charge in each time constant bin.

    # Up-the-ramp exposure; Repeat N times.
    for i in range(0,N):
        Qij,Qi,Ei = persistence_const(x,rho,E0[i],Qij)
        E.append(Ei)
        Q.append(Qi)

    # Dark
    N = 2000
    T = 5.0 #sec
    x = T/tau

    tframe = np.concatenate([tframe,np.array(range(0,N))*T+tframe[-1]])
    Ed0 = np.zeros(N)
    E0  = np.concatenate([E0,Ed0])
    for i in range(0,N):
        Qij,Qi,Ei = persistence_const(x,rho,Ed0[i],Qij)
        E.append(Ei)
        Q.append(Qi)

    fs  = 18
    fig = plt.figure(figsize=(10,7))
    ax  = fig.add_subplot(212)
    ax.plot(tframe,Q,".")
    plt.ylabel("trapped signal [e-]",fontsize=fs)
    plt.xlabel("time [sec]",fontsize=fs)
    plt.ylim(0,400.0)
    plt.xlim(0,35000.0)
    plt.tick_params(labelsize=fs)

    ax3 = fig.add_subplot(211)
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
