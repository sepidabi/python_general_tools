from My_functions import read_polcube as r
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as ply
import pickle
from PIL import Image

dir='/srv/scratch/crobu/sst/2017.04.20/'
#dir='/scratch/crobu/2017.04.20/'
f0=dir+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected.fcube'
f1=dir+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected_sp.fcube'



ncl=7
cube=r(f0,f1)
nt=cube[:,0,0:21,:,:].shape[0]
nw=cube[:,0,0:21,:,:].shape[1]
ny=cube[:,0,0:21,:,:].shape[2]
nx=cube[:,0,0:21,:,:].shape[3]

dat1=np.zeros((nw,ny,nx))
dat_profs=np.zeros((nw,7,1))


im_frame = Image.open(dir+'mask_chromis.png')
mask= np.array(im_frame.getdata()).reshape(ny,nx)



# for x in range(nx):
#     print(x,'/',nx-1)
#     for y in range(ny):
#         if (mask[y,x] == 0):
#             dat1[:,y,x]=cube[26,0,0:21,y,x]/cube[26,0,0:21,y,x].max()

dat1=cube[26,0,0:21,:,:]      

# for x in range(nx/6):
#     for y in range(ny):
#         dat_6prof[:,y,x]=cube[26,0,0:21,380,493]
#         dat_6prof[:,y,x+nx/6]=cube[26,0,0:21,472,1472]
#         dat_6prof[:,y,x+nx/6*2]=cube[26,0,0:21,89,1638]
#         dat_6prof[:,y,x+nx/6*3]=cube[26,0,0:21,518,1053]
#         dat_6prof[:,y,x+nx/6*4]=cube[26,0,0:21,114,1449]
#         dat_6prof[:,y,x+nx/6*5]=cube[26,0,0:21,637,332]

dat_profs[:,0,0]=cube[26,0,0:21,360,799]/cube[26,0,0:21,380,493].max()
dat_profs[:,1,0]=cube[26,0,0:21,354,1248]/cube[26,0,0:21,472,1472].max()
dat_profs[:,2,0]=cube[26,0,0:21,89,1638]/cube[26,0,0:21,89,1638].max()
dat_profs[:,3,0]=cube[26,0,0:21,518,1053]/cube[26,0,0:21,518,1053].max()
dat_profs[:,4,0]=cube[26,0,0:21,114,1499]/cube[26,0,0:21,114,1499].max()
dat_profs[:,5,0]=cube[26,0,0:21,637,332]/cube[26,0,0:21,637,332].max()
dat_profs[:,6,0]= np.zeros(nw)

nnx=1
nny=7
mm= np.zeros((nny*nnx))


dat =dat_profs.reshape(nw,nnx*nny)
clusterer = KMeans(n_clusters=ncl).fit(np.transpose(dat))

dat2 = dat1.reshape(nw,nx*ny)
dat2_labels = clusterer.predict(np.transpose(dat2))
x=clusterer.cluster_centers_

# dat =dat1.reshape(nw,nx*ny)
# mm[:]=KMeans(n_clusters=ncl).fit(np.transpose(dat)).labels_
# x=KMeans(n_clusters=ncl).fit(np.transpose(dat)).cluster_centers_





f = open(dir+'clusters.pckl', 'wb')
pickle.dump(dat2_labels, f)
f.close()

f = open(dir+'labels.pckl', 'wb')
pickle.dump(x, f)
f.close()


# f = open(dir+'clusters.pckl', 'rb')
# dat2_labels = pickle.load(f)
# f.close()

# f = open(dir+'labels.pckl', 'rb')
# x = pickle.load(f)
# f.close()

    
ply.imshow(dat2_labels.reshape(ny,nx), origin='lower')

ply.figure()

for i in range(ncl):
    #ply.plot(x[i,:],(i*2/10.,i*3/10.,i*5/10.),label = str(i))
    ply.plot(x[i,:],label = str(i))
ply.legend()

ply.show()
