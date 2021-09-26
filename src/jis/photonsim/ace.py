import math
import numpy as np
import pyfftw
import sys
import matplotlib.pyplot as plt

def calc_ace(rg, N, T, ace):
    """
    This function makes a one-dimensional attitude control error (ACE) data.
    The ACE is assumed to have a PSW which consists of a power-law function
    and some Lorenzian-type peaks (PS disturbance). The returned ace data (acedata)
    is normalized with its standard deviation. The generated PSD is also returned
    as 'psdn', but this is not normalized as 'acedata'.

    Args:
        rg  (numpy.random.Generator): Random generator.
        N   (int)  : Number of time grids.
        T   (float): Range of time to simulate.
        ace (dict) : Parameterset about ace loaded from the ace json file.

    Returns:
        acedata (ndarray): Time-series data of the calculated ace normalized with stddev.
        psdn    (ndarray): Calculated PSD data (Not normalized; just for checking).

    """

    dT = T/N # Time step.

    # Set parameters and PSD function
    title = ace['title']
    # print("used model {0},{1},{2}".format(title,ace['model'],ace['note']))
    if title != 'ace001':
        print("Not supported title {0}".format(title))
        sys.exit()
    else:
        f0   = ace['f0']  # Reference frequency to define the power-law PSD.
        alp  = ace['alp'] # Factor which determines the fraction of the gaussian component.
        bet  = ace['bet']
        nmax = ace['n']   # Number of the PS disturbance peaks (see below).

        # Settings for the PS disturbance (peaks of some specific frequencies).
        if nmax > 0:
            P = np.empty(nmax) # Power
            F = np.empty(nmax) # Frequency
            H = np.empty(nmax) # Half-width half-maximum
            for i in range(nmax):
                pi = "p{}".format(i+1)
                hi = "h{}".format(i+1)
                fi = "f{}".format(i+1)
                P[i] = ace[pi]
                H[i] = ace[hi]
                F[i] = ace[fi]

        """
        -----------------------------------------------
        Set Fourier domain ACE

        ACE(t) should be a real (not complex) function.
        Fourier domain ACE F(f) should...
          F(0 or N/2 ) should be real
          F(N-f) = conjugete( F(f) )
            then
              ReF(-f) = ReF(f) ,  ImF(-f) = -ImF(f)
        -----------------------------------------------
        """

        # initialize data for fft
        data = pyfftw.zeros_aligned((N),dtype='complex128')

        # independent data (i = 0,...,N/2) ######################
        t   = np.arange(0, int(N/2)+1, dtype='int')

        f   = t/N/dT                 # Freq. grid.
        psd = 1/f0/(1+(f/f0)**2.)    # PSD of the continuum comp.

        # Adding PS-disturbance peaks (Lorenzian function).
        if nmax > 0:
            for i in range(nmax):
                psd = psd + P[i]/H[i]/(1+((f-F[i])/H[i])**2)

        s0    = np.sqrt(psd)         # Amplitude without noise.
        sig   = s0**bet              # Sigma of the gaussian noise in the next line.
        noise = rg.normal(scale=sig) # Gaussian noise.
        s     = s0 + alp * noise     # Amplitude with noise.

        # Making phase th (=theta).
        th = rg.uniform(size=np.size(s))*2.*np.pi

        ## To make the data real values, set some phases to be zero.
        th[0] = 0.
        if N%2 == 0:
            th[int(N/2)] = 0.

        # Making Fourier components.
        data.real[t] = s*np.cos(th)
        data.imag[t] = s*np.sin(th)

        # dependent data (i = N/2+1,...,N-1) ####################
        t = np.arange(0, int(N/2)-1, dtype='int')
        tt = N - t - 1
        data.real[tt] =  data.real[t+1]
        data.imag[tt] = -data.imag[t+1]

        # iFFT
        ft = pyfftw.interfaces.numpy_fft.ifft(data)

        # statistics
        mean    = np.mean( ft.real )
        std     = np.std( ft.real, ddof=1) # sqrt( SUM((x-<x>)**2)/(N-1) )

        # Making output data.
        acedata = ( ft.real - mean ) / std # Centering and normalizing with std.
        psdn    = data.real*data.real+data.imag*data.imag
                                           # PSD of the result (not normalized).

    return acedata, psdn


def calc_dummy_ace(rg, N, T, ace):
    """
    This function generates a one-dimensional dummy attitue control error (ACE) data.
    The arguents are the same as `calc_ace` but not used except for the number of elements, `N`.
    The generated ACE is just a zero-padded one-dimensional array.

    Args:
        rg  (numpy.random.Generator): Random generator.
        N   (int)  : Number of time grids.
        T   (float): Range of time to simulate.
        ace (dict) : Parameter set about ace loaded from the ace json file.

    Returns:
        dummyace (ndarray): A dummy ACE, zero-padded array with the size of `N`.

    """
    return np.zeros(N)


def plot_ace(N, T, data, psdn, plotfile):
    """
    This function makes a plot of the ace data.
    The plot is saved as a file, 'plotfile'.

    Args:
        N    (int)    : Number of time grids.
        T    (float)  : Range of time simulated.
        data (ndarray): Simulated ACE data.
        psdn (ndarray): Simulated PSD data.
        plotfile (str): Output file name.

    Returns:
        None.

    """

    dT  = T/N # time step
    Ta  = np.arange(0, N*dT, dT)
    sig = np.std(data, ddof=1)

    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    p1  = ax1.plot(Ta, data)

    ax1.set_ylabel('Atitude Control Error [arcsec]')
    ax1.set_xlabel('Time [sec]')
    ax1.text(dT*20,3.0,"Sigma={:.2f}".format(sig))
    ax1.set_ylim(bottom=-4,top=4)

    FR = np.arange(0, int(N/2+1))
    FR = FR / (N*dT)

    ax2 = fig.add_subplot(2,1,2)
    p2  = ax2.plot(FR,psdn[0:int(N/2+1)], marker='o',\
                   markersize=1,linestyle='None')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel('PSD')
    ax2.set_xlabel('Frequency [Hz]')

    plt.subplots_adjust(hspace=0.4)

    fig.savefig(plotfile,transparent=True)

    plt.show()


def ace2d(x_scale, y_scale, pix_scale, N, xdata, ydata, xN):
    """
    This function returns 2D ace data array with considering
    the ace data scale and the pixel scale.

    Args:
        x_scale   (float)  : Scale of xdata.
        y_scale   (float)  : Scale of ydata.
        pix_scale (float)  : Pixel scale of the output image.
        N         (int)    : Array size of the output image (N x N).
        xdata     (ndarray): X-axis simulated ACE array.
        ydata     (ndarray): Y-axis simulated ACE array
        xN        (int)    : Simulated ACE array size.

    Returns:
        data      (ndarray): 2D array of the ace data.

    """

    data = np.zeros( (N,N) , dtype=int)

    for i in range(xN):
        x  = xdata[i]*x_scale/pix_scale
        y  = ydata[i]*y_scale/pix_scale
        xi = math.floor(x + 0.5 + N/2 )
        yi = math.floor(y + 0.5 + N/2 )
        data[yi,xi] = data[yi,xi] + 1

    return data
