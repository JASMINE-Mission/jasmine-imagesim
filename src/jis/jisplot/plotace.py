"""plotting for ACE

"""

import matplotlib.pyplot as plt

def trajectory(xdata,ydata,save=False,fs=16):
    """
    Summary
    --------
    plotting for trajectory

    Parameters
    --------------
    xdata : ndarray
            trajectory for x axis
    ydata : ndarray
            trajectory for y axis
    
    Returns
    -------------
    
    """

    fig=plt.figure(figsize=(7,7))
    ax=fig.add_subplot(111)
    ax.plot(xdata,ydata,".",alpha=0.3)
#    ax.set_aspect(1.0/ax.get_data_ratio())
    plt.xlabel("Pixel scale",fontsize=fs)
    plt.ylabel("Pixel scale",fontsize=fs)
    plt.tick_params(labelsize=fs)
    if save:
        plt.savefig(save, bbox_inches="tight", pad_inches=0.0)
    else:
        plt.show()
        
    return 
