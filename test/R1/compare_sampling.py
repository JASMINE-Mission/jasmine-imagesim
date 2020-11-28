import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
with h5py.File("pixcube_samp1.h5", "r") as f:
    data1=f["data/pixcube"].value
with h5py.File("pixcube_samp2.h5", "r") as f:
    data2=f["data/pixcube"].value
with h5py.File("pixcube_samp3.h5", "r") as f:
    data3=f["data/pixcube"].value
with h5py.File("pixcube_samp4.h5", "r") as f:
    data4=f["data/pixcube"].value

nx,ny,nt=np.shape(data4)
gg=np.ones(nx)
gx=np.array(range(0,nx))
X=np.array([gg]).T*gx
Y=X.T
#data4=np.ones(np.shape(data4))

centX1=np.sum(X*data1[:,:,0])/np.sum(data1[:,:,0])
centY1=np.sum(Y*data1[:,:,0])/np.sum(data1[:,:,0])

centX2=np.sum(X*data2[:,:,0])/np.sum(data2[:,:,0])
centY2=np.sum(Y*data2[:,:,0])/np.sum(data2[:,:,0])

centX3=np.sum(X*data3[:,:,0])/np.sum(data3[:,:,0])
centY3=np.sum(Y*data3[:,:,0])/np.sum(data3[:,:,0])

centX4=np.sum(X*data4[:,:,0])/np.sum(data4[:,:,0])
centY4=np.sum(Y*data4[:,:,0])/np.sum(data4[:,:,0])

print(np.sqrt((centX4-centX3)**2+(centY4-centY3)**2))
print(np.sqrt((centX4-centX2)**2+(centY4-centY2)**2))
print(np.sqrt((centX4-centX1)**2+(centY4-centY1)**2))


print(np.shape(data1))
a=plt.imshow((data4[:,:,0]-data1[:,:,0])/data4[:,:,0])
plt.colorbar(a)
plt.savefig("c41.png")
plt.show()
