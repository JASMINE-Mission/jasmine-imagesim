import h5py
import matplotlib.pyplot as plt 
import numpy as np
with h5py.File("image.h5","r") as f:
    cube=f["data/pixcube"].value
print(np.shape(cube))
plt.imshow(cube[0:200,0:200,0])
plt.show()
s=50
e=100
lc=np.sum(cube[s:e,s:e,:],axis=(0,1))
norm=np.median(lc)
lc=lc/norm
plt.plot(lc)
plt.show()
