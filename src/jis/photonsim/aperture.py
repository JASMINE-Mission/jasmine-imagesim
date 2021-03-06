import numpy as np
import math

def calc_aperture(N,EPD,Robs,Tsp):
    """
    This function calculates the aperture pattern and its area.

    The aperture pattern consists of a circular aperture
    with an entrance pupil diameter of EPD, a central 
    circular obscuration with a radius of Robs, and
    a 3-vane spider pattern with a thickness of Tsp. 

    The output pattern consists of N x N apt-cells with a value of 0 or 1.
    The parameters (EPD, Robs, and Tsp) are in "apt-cell".
    But in the current simulation, we are assuming a scale of 1 mm/apt-cell.

    Args:
        N    (int)  : Number of apt-cells of the output aperture pattern (apt-cell).
        EPD  (float): Entrance pupil diameter (apt-cell).
        Robs (float): Radius of the central obscuration (apt-cell).
        Tsp  (float): Thickness of the spider (apt-cell).
        *** Scale is assumed to be 1 mm/apt-cell!!! ***

    Returns:
        data (ndarray): N x N data array of the aperture pattern.
        S    (float)  : Total area (number of apt-cells=mm^2) of unmasked region.

    Examples:
        from jis.photonsim.aperture import calc_aperture

        N    = 100 # Number of apt-cells of the output aperture pattern.
        EPD  = 80  # Entrance pupil diameter in apt-cell.
        Robs = 20  # Central obscuration radius in apt-cell.
        Tsp  = 3   # Spider thickness in apt-cell.

        data, S = calc_aperture(N, EPD, Robs, Tsp)


    """
    # initialize data
    data=np.zeros((N,N),dtype=np.float64) # 

    # Set aperture
    a = math.sqrt(3)/2
    S=0
    for i in range(N):
        y = i - N/2 # y axis,  origin at N/2
        for j in range(N):
            x   =    j - N/2  # x axis
            rsq =  x*x + y*y  # square of distance from the center 
            x1  = -x/2 + y*a  # x1,y1 is rotated -120 degree 
            y1  = -x*a - y/2
            x2  = -x/2 - y*a  # x2,y2 is rotated 120 degree 
            y2  =  x*a - y/2
            apt = 0.0 
            if rsq <= (EPD/2.)**2  and rsq > Robs**2 : 
                apt = 1.0 
                if x  >0 and y  > -Tsp/2 and y  < Tsp/2: # in a spider
                    apt = 0.0
                if x1 >0 and y1 > -Tsp/2 and y1 < Tsp/2: # in a spider
                    apt = 0.0
                if x2 >0 and y2 > -Tsp/2 and y2 < Tsp/2: # in a spider
                    apt = 0.0
            #NOTE: if Rob==0 and Tsp==0 then apt=1 for r==0 and y==0 
            data[i,j] = apt 
            S = S + apt
    return data, S
