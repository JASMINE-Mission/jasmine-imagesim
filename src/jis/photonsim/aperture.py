import numpy as np
import math
def calc_aperture(N,EPD,Robs,Tsp):
    # initialize data
    data=np.zeros((N,N),dtype=np.float64) # 
    # Set aperture
    a = math.sqrt(3)/2
    S=0
    for i in range(N):
        y = i - N/2 # y axis,  origin at N/2
        for j in range(N):
            x = j - N/2   # x axis
            r = x*x + y*y # square of distance from the center 
            x1 = -x/2 + y*a  # x1,y1 is rotated -120 degree 
            y1 = -x*a - y/2
            x2 = -x/2 - y*a  # x2,y2 is rotated 120 degree 
            y2 =  x*a - y/2
            apt = 0.0 
            if r <= (EPD/2.)**2  and r > Robs**2 : 
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
