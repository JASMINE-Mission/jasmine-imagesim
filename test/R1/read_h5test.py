import h5py
import matplotlib.pyplot as plt 
import numpy as np
with h5py.File("image.h5","r") as f:
    cube=f["data/pixcube"].value
print(np.shape(cube))
plt.plot(np.sum(cube[0:,0:,:],axis=(0,1)))
plt.show()
