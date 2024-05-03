'''
==============
3D quiver plot
==============

Demonstrates plotting directional arrows at points on a 3d meshgrid.
'''

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from scipy.io.idl import readsav 

vect = readsav('/srv/scratch/crobu/sst/2013.07.15/b_field_rotated_coupled_allcolumn.csav')

step = 10 
taup=30

az=vect.azim_local[taup,:,:]


inc=(np.arctan(vect.bz / np.sqrt(vect.bx**2+vect.by**2)))[30,0:100,0:100]
inc= [inc[i] for i in xrange(0, len(inc), 5)]
B=(np.sqrt(vect.bx**2+vect.by**2+vect.bz**2))[30,0:100,0:100]
B= [B[i] for i in xrange(0, len(B), 5)]

fig = plt.figure()
#ax = fig.gca(projection='3d')
ax = fig.add_subplot(111, projection='3d')

# Make the grid
#x = np.arange(0, 10, 1)
#y = np.arange(0, 10, 1)
#z = (np.arange(0, 10, 1))*0


x, y, z = np.meshgrid(np.arange(0, 100, 5),
                      np.arange(0, 100, 5),
                      np.arange(0, 100, 5)*0)
                     

# Make the direction data for the arrows

u = B*np.cos(az)*np.cos(inc)
v = B*np.sin(az)*np.cos(inc)
w = B*np.sin(inc)


ax.quiver(x, y, z, u, v, w)
#ax.set_xlim3d(np.min(u),np.max(u))
#ax.set_ylim3d(np.min(v),np.max(v))
#ax.set_zlim3d(np.min(w),np.max(w))
ax.set_xlim([0,100])
ax.set_ylim([0,100])
ax.set_zlim([-1.5,0])

plt.show()
