import numpy as np
from My_functions import *
import matplotlib.pyplot as plt
import scipy.io as sc
from PIL import Image, ImageDraw


def plot_prof(xc,yc,side):
    x0=xc-side
    x1=xc+side
    y0=yc-side
    y1=yc+side
    ax2.plot(np.median(np.median(cube_cak[:,y0:y1,x0:x1], axis=1), axis=1), label=str(side))

def plot_rectangle(xc,yc,side,col):
    x0=xc-side
    x1=xc+side
    y0=yc-side
    y1=yc+side
    #draw.rectangle([(x0, y0), (x1, y1)],outline=col)
    ax.plot



dir='/srv/scratch/crobu/sst/2017.04.20/'
file0=dir+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected.fcube'
file1=dir+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected_sp.fcube'

tt=26

cube_cak=read_polcube(file0,file1)[tt,0,0:21,:,:]

file2=dir+'crispex.stokes.8542.09:40:25.time_corrected_CHROMIS.fcube'
file3=dir+'crispex.stokes.8542.09:40:25.time_corrected_CHROMIS_sp.fcube'
cube_8542=read_polcube(file2,file3)



spect_pos_k=sc.readsav(dir+'spectfile.3950.idlsave',verbose=True,python_dict=True)['spect_pos'][0:21]

nt= cube_8542.shape[0]
nw= spect_pos_k.shape[0]
nw_8542=cube_8542.shape[2]
ny= cube_8542.shape[3]
nx= cube_8542.shape[4]


deriv= np.zeros((nw,ny,nx))
deriv_sec= np.zeros((nw,ny,nx))
number_dzero= np.zeros((ny,nx))
single_loc=  np.zeros((nw,ny,nx))
v= np.zeros((nw_8542,ny,nx))
im_frame = Image.open('mask_chromis.png')
mask= np.array(im_frame.getdata()).reshape(ny,nx)


# for x in range(nx):
#     print(x,'/',nx-1)
#     for y in range(ny):
#         deriv[:,y,x]= derivative(spect_pos_k, cube_cak[:,y,x])
#         deriv_sec[:,y,x]= derivative(spect_pos_k, deriv[:,y,x])
#         dum=np.where(deriv[:,y,x] == 0.)
#         a=dum[0]
#         number_dzero[y,x]=len(a)
#         if (mask[y,x] == 0):
#             v[:,y,x]= np.sum(cube_8542[:,3,:,y,x], axis=0) #cube_8542[tt,3,:,y,x]
#             if (number_dzero[y,x] == 3) and (deriv_sec[a[0],y,x] > 0.) and  (deriv_sec[a[1],y,x] < 0.):
#                 single_loc[:,y,x] = cube_cak[:,y,x]

xc=1217
yc=385
side=40

x0=xc-side
x1=xc+side

y0=yc-side
y1=yc+side

fig, (ax1, ax2) = plt.subplots(1, 1, sharey=True)
ax1.imshow(cube_cak[10,:,:])


# for i in range(5):
#     col= 50*i
#     a=plot_rectangle(xc,yc,10*(i+1),col)


for i in range(5):
    a=plot_prof(xc,yc,10*(i+1))
    
plt.legend()

fig.savefig('prova.eps', format='eps', dpi=1000)




# plt.plot(np.median(np.median(cube_cak[:,y0:y1,x0:x1], axis=1), axis=1))
# plt.figure()
# plt.imshow(cube_cak[10,:,:], origin='lower', cmap='gist_gray')


#plt.imshow(single_loc[0,:,:], cmap='gist_gray', origin='lower')
#plt.contour(v[6,:,:])

