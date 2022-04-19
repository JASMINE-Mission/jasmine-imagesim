import astropy.io.ascii as asc
import matplotlib.pylab as plt
def plot_variability(variability,filenames_starplate,tday):
    """plot  variability
    
    Args:
        variability: variability object
        filenames_starplate: filenames list of starplate
        tday: time bin in the unit of day

    """
    for line in asc.read(filenames_starplate):
        varsw, injlc, b=variability.read_var(tday,line['star index'])
        if varsw:
            plt.plot(tday,injlc)
            plt.savefig("variability_input"+"_"+str(line['star index'])+".png")
            plt.clf()
