## I dont impose any external cluster_center


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

a=[1120,1500]
b=[300,487]

time_step=26

# a=[900,1500]
# b=[300,600]

ncl=4
cube=r(f0,f1)
nt=cube[:,0,0:21,b[0]:b[1],a[0]:a[1]].shape[0]
nw=cube[:,0,0:21,b[0]:b[1],a[0]:a[1]].shape[1]
ny=cube[:,0,0:21,b[0]:b[1],a[0]:a[1]].shape[2]
nx=cube[:,0,0:21,b[0]:b[1],a[0]:a[1]].shape[3]



dat1=np.zeros((nw,ny,nx))

# im_frame = Image.open(dir+'mask_chromis.png')
# mask= np.array(im_frame.getdata()).reshape(ny,nx)



# for x in range(nx):
#     print(x,'/',nx-1)
#     for y in range(ny):
#         if (mask[y,x] == 0):
#             dat1[:,y,x]=cube[26,0,0:21,y,x]#/cube[26,0,0:21,y,x].max()


dat1=cube[time_step,0,0:21,b[0]:b[1],a[0]:a[1]]


mm= np.zeros((ny*nx))


dat =dat1.reshape(nw,nx*ny)
cluster=KMeans(n_clusters=ncl).fit(np.transpose(dat))
mm[:]=cluster.labels_
x=cluster.cluster_centers_





f = open(dir+'clusters.pckl', 'wb')
pickle.dump(mm, f)
f.close()

f = open(dir+'labels.pckl', 'wb')
pickle.dump(x, f)
f.close()


# f = open(dir+'clusters.pckl', 'rb')
# mm = pickle.load(f)
# f.close()

# f = open(dir+'labels.pckl', 'rb')
# x = pickle.load(f)
# f.close()

colors = ply.cm.jet(np.linspace(0,1,ncl))

ply.imshow(mm.reshape(ny,nx), origin='lower', cmap='jet')

ply.figure()

for i in range(ncl):
    ply.plot(x[i,:],label = str(i),color=colors[i])
ply.legend()

# ply.figure()
# ply.imshow(cube[time_step,0,10,b[0]:b[1],a[0]:a[1]], origin='lower', cmap='gist_gray')

ply.show()
