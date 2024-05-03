from My_functions import *
from My_functions import read_polcube as r
from My_functions import read_lp_cube as read_lp_cube
import numpy as np
from sklearn.cluster import KMeans
import matplotlib 
import matplotlib.pyplot as plt
import pickle
from PIL import Image
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors
import pyana as pyana
import sys
import os
from astropy.io import fits
from My_functions import read_polcube as r
from matplotlib import rc
import pylab as pl
import matplotlib.ticker as ticker

def mfits(file):
    a = fits.open(file, 'readonly')
    return np.float32(a[0].data[:])



dir='/srv/scratch/crobu/sst/2017.04.20/'
f0=dir+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected.fcube'
f1=dir+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected_sp.fcube'

a=[1120,1500]
b=[300,487]

time_step=26

ncl=4
cube=r(f0,f1)
nt=cube[:,0,0:21,b[0]:b[1],a[0]:a[1]].shape[0]
nw=cube[:,0,0:21,b[0]:b[1],a[0]:a[1]].shape[1]
ny=cube[:,0,0:21,b[0]:b[1],a[0]:a[1]].shape[2]
nx=cube[:,0,0:21,b[0]:b[1],a[0]:a[1]].shape[3]

# dat1=np.zeros((nw,ny,nx))
# dat1=cube[time_step,0,0:21,b[0]:b[1],a[0]:a[1]]

# dat =dat1.reshape(nw,nx*ny)
# clusterer = KMeans(n_clusters=ncl).fit(np.transpose(dat))

# dat2_labels=np.zeros((nx*ny,nt))
# x= np.zeros((ncl,nw,nt))


# for tt in range(nt):
#     dat2 =cube[tt,0,0:21,b[0]:b[1],a[0]:a[1]]
#     dat2 =dat2.reshape(nw,nx*ny)
#     dat2_labels[:,tt] = clusterer.predict(np.transpose(dat2))
#     x[:,:,tt]=clusterer.cluster_centers_

# imm=dat2_labels.reshape(ny,nx,nt)

# file = open(dir+'dat2_labels.pckl', 'wb')
# pickle.dump(dat2_labels, file)
# file.close()

file = open(dir+'dat2_labels.pckl', 'rb')
dat2_labels = pickle.load(file)
file.close()


#ply.imshow(dat2_labels.reshape(ny,nx,nt)[:,:,26], origin='lower')


rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching  helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]


mu=0.33
cont = {'cak':2.31863451e-05/get_limbdarkening(mu,'3950')}


calib_3934 =  mfits(dir+'calib.3934.fits')
calib_3968 =  mfits(dir+'calib.3968.fits')
wl=((pyana.fzread(dir+'wav.3950.f0'))['data'])*1e10 #in angstrom
w0=3933.6640
wlc=((wl[0:21]+ calib_3934[1])*calib_3934[2])
vdopp= (1-w0/wlc)*299792.458

params = {'figure.figsize': [3.46457,3.46457], 'font.size':8, 'axes.labelsize': 8,
          'font.family':matplotlib.font_manager.FontProperties(fname='/home/crobu/HELR45W.ttf').get_name(),
      }
matplotlib.rcParams.update(params)


f = plt.figure()
f, ax = plt.subplots(3, 2)


ax1=ax[0,0]
ax2=ax[0,1]
ax3=ax[1,0]
ax4=ax[1,1]
ax5=ax[2,0]
ax6=ax[2,1]

colors = ('c','limegreen','r','k')
level=[0,1,2,2.9]
steps=[6,41,70]

assex=np.arange(a[0],a[1],1)*0.038
assey=np.arange(b[0],b[1],1)*0.038

i=0
for ax in [ax2,ax4,ax6]:
    ax.pcolormesh(assex,assey,cube[steps[i],0,10,b[0]:b[1],a[0]:a[1]]**(-1), cmap='gist_gray',shading='gouraud')
    ax.axis('image')
    ax.tick_params(axis='both', direction='out',length=2)
    ax.set_yticklabels([])
    if (ax==ax2 or ax==ax4):
        ax.set_xticklabels([])
    if (ax==ax6):
        ax.set_xlabel(r'$\mathrm{X[arcsec]}$')
    ax.contour(np.arange(a[0],a[1],1)*0.038,np.arange(b[0],b[1],1)*0.038,dat2_labels.reshape(ny,nx,nt)[:,:,steps[i]],levels=[2.9], colors=colors[1], linewidths=1.5)
    ax.contour(np.arange(a[0],a[1],1)*0.038,np.arange(b[0],b[1],1)*0.038,dat2_labels.reshape(ny,nx,nt)[:,:,steps[i]],levels=[2],  colors=colors[0], linewidths=1.)
    ax.contour(np.arange(a[0],a[1],1)*0.038,np.arange(b[0],b[1],1)*0.038,dat2_labels.reshape(ny,nx,nt)[:,:,steps[i]],levels=[1], colors=colors[2], linewidths=1.)
    ax.contour(np.arange(a[0],a[1],1)*0.038,np.arange(b[0],b[1],1)*0.038,dat2_labels.reshape(ny,nx,nt)[:,:,steps[i]],levels=[0], colors=colors[3], linewidths=1.)
    i +=1

vdopp=read_lp_cube(dir+'vdopp_cak3.fcube')



i=0
for ax in [ax1,ax3,ax5]:
    ax.pcolormesh(assex,assey,vdopp[steps[i],b[0]:b[1],a[0]:a[1]]*(-1)+4.3,cmap='RdBu', vmin=-20, vmax=20,shading='gouraud',clip_on=True)
    ax.axis('image')
    ax.tick_params(axis='both', direction='out',length=2)
    if (ax==ax1 or ax==ax3):
        ax.set_xticklabels([])
    if (ax==ax5):
        ax.set_xlabel(r'$\mathrm{X[arcsec]}$')
    ax.set_ylabel(r'$\mathrm{Y[arcsec]}$')
    i +=1

f.set_tight_layout(True)



pp = PdfPages('clustering_time_series.pdf')
pp.savefig(bbox_inches='tight', pad_inches=0.05)
pp.close()

