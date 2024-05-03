from My_functions import read_polcube as r
import numpy as np
import matplotlib.pyplot as ply
import pickle
from PIL import Image
from numpy import array
from scipy.cluster.vq import vq


dir='/srv/scratch/crobu/sst/2017.04.20/'
#dir='/scratch/crobu/2017.04.20/'
f0=dir+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected.fcube'
f1=dir+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected_sp.fcube'

cube=r(f0,f1)
nt=cube[:,0,0:21,:,:].shape[0]
nw=cube[:,0,0:21,:,:].shape[1]
ny=cube[:,0,0:21,:,:].shape[2]
nx=cube[:,0,0:21,:,:].shape[3]

mm= np.zeros((ny*nx))
dat1=np.zeros((nw,ny,nx))

im_frame = Image.open(dir+'mask_chromis.png')
mask= np.array(im_frame.getdata()).reshape(ny,nx)

dat1=cube[26,0,0:21,:,:]

# for x in range(nx):
#     print(x,'/',nx-1)
#     for y in range(ny):
#         #if (mask[y,x] == 0):
#             dat1[:,y,x]=cube[26,0,0:21,y,x]
            
dat =dat1.reshape(nw,nx*ny)
#coords=[(380,493),(637,332)]
coords=[dat1[:,380,493],dat1[:,637,332]]

x=vq(dat1,coords,)
