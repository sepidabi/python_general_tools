from My_functions import *
from My_functions import read_polcube as r
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib 
import pickle
from PIL import Image
from matplotlib.backends.backend_pdf import PdfPages
import pyana as pyana
import sys
import os
from astropy.io import fits
from My_functions import read_polcube as r
from matplotlib import rc
import pylab as pl
import matplotlib.gridspec as gridspec
def mfits(file):
    a = fits.open(file, 'readonly')
    return np.float32(a[0].data[:])

 
rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching  helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]  



a=[1120,1500]
b=[300,487]
time_step=26

mu=0.33
cont = {'cak':2.31863451e-05/get_limbdarkening(mu,'3950')}


dir='/srv/scratch/crobu/sst/2017.04.20/'
f0=dir+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected.fcube'
f1=dir+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected_sp.fcube'
ncl=4
cube=r(f0,f1)
nt=cube[:,0,0:21,b[0]:b[1],a[0]:a[1]].shape[0]
nw=cube[:,0,0:21,b[0]:b[1],a[0]:a[1]].shape[1]
ny=cube[:,0,0:21,b[0]:b[1],a[0]:a[1]].shape[2]
nx=cube[:,0,0:21,b[0]:b[1],a[0]:a[1]].shape[3]
dat1=np.zeros((nw,ny,nx))
dat1=cube[time_step,0,0:21,b[0]:b[1],a[0]:a[1]]


dir='/srv/scratch/crobu/sst/2017.04.20/'
f = open(dir+'clusters.pckl', 'rb')
dat2_labels = pickle.load(f)
f.close()

f = open(dir+'labels.pckl', 'rb')
x = pickle.load(f)
f.close()

calib_3934 =  mfits('calib.3934.fits')
calib_3968 =  mfits('calib.3968.fits')
wl=((pyana.fzread(dir+'wav.3950.f0'))['data'])*1e10 #in angstrom
w0=3933.6640
wlc=((wl[0:21]+ calib_3934[1])*calib_3934[2])
vdopp= (1-w0/wlc)*299792.458-4.3



params = {'figure.figsize': [3.46457,3.46457*1.3], 'font.size':8, 'axes.labelsize': 8,
          'font.family':matplotlib.font_manager.FontProperties(fname='/home/crobu/HELR45W.ttf').get_name(),
      }
matplotlib.rcParams.update(params)


# Rig=7
# Col=4

f = plt.figure()
# ax1 = plt.subplot2grid((Rig,Col), (0, 0), colspan=4, rowspan=3)
# ax2 = plt.subplot2grid((Rig,Col), (3, 0), colspan=2, rowspan=2)
# ax3 = plt.subplot2grid((Rig,Col), (3, 2), colspan=2, rowspan=2)
# ax4 = plt.subplot2grid((Rig,Col), (5, 0), colspan=2, rowspan=2)
# ax5 = plt.subplot2grid((Rig,Col), (5, 2), colspan=2, rowspan=2)


outer =gridspec.GridSpec(2,1)

AXA = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec = outer[0],hspace = 0.5)
ax1  = plt.subplot(AXA[0,:])

AXB = gridspec.GridSpecFromSubplotSpec(2,2,subplot_spec = outer[1])
ax2 = plt.subplot(AXB[0,0])
ax3 = plt.subplot(AXB[0,1])
ax4 = plt.subplot(AXB[1,0])
ax5 = plt.subplot(AXB[1,1])



colors = ('c','limegreen','r','k')
level=[0,1,2,2.9]

ax1.pcolormesh(np.arange(a[0]-100,a[1]-100,1)*0.038,np.arange(b[0],b[1],1)*0.038, dat1[10,:,:]**(-1),  cmap='gist_gray',shading='gouraud')
ax1.axis('image')
ax1.tick_params(axis='both', direction='out',length=2)
#ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
#ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
ax1.set_xlabel(r'$\mathrm{X[arcsec]}$')
ax1.set_ylabel(r'$\mathrm{Y[arcsec]}$')
ax1.contour(np.arange(a[0]-100,a[1]-100,1)*0.038,np.arange(b[0],b[1],1)*0.038,dat2_labels.reshape(ny,nx),levels=[2.9], colors=colors[3], linewidths=3.5)
ax1.contour(np.arange(a[0]-100,a[1]-100,1)*0.038,np.arange(b[0],b[1],1)*0.038,dat2_labels.reshape(ny,nx),levels=[2],  colors=colors[2], linewidths=1.)
ax1.contour(np.arange(a[0]-100,a[1]-100,1)*0.038,np.arange(b[0],b[1],1)*0.038,dat2_labels.reshape(ny,nx),levels=[1], colors=colors[1], linewidths=1.4)
ax1.contour(np.arange(a[0]-100,a[1]-100,1)*0.038,np.arange(b[0],b[1],1)*0.038,dat2_labels.reshape(ny,nx),levels=[0], colors=colors[0], linewidths=1.4)
x= x/(1e-7)
i = 0
for ax in [ax2,ax3,ax4,ax5]:
    ax.plot(vdopp,x[i,:]/calib_3934[0]*get_limbdarkening(mu,'3950'),label = str(i),color=colors[i],linewidth=2)
    ax.set_xticks([-50,0,50])
    ax.set_xticks(np.arange(-50,50,10), minor=True)
    if (ax==ax4 or ax==ax5):
        ax.set_xlabel(r'$\mathrm{V}_{\mathrm{los}}[\mathrm{km} \, \mathrm{s}^{-1}]$', labelpad=1)
    if (ax==ax2 or ax==ax4):
        ax.set_ylabel(r'$\mathrm{Intensity}$') # (\mathrm{erg}\, \mathrm{s}^{-1} \, \mathrm{cm}^{-2} \, \mathrm{Hz}^{-1} \, \mathrm{ster}^{-1})$')
    if (ax==ax2 or ax==ax3):
        ax.set_xticklabels([])
    i += 1

#AX.set_tight_layout(True)


pp = PdfPages('clustering.pdf')
pp.savefig(bbox_inches='tight', pad_inches=0.03)
pp.close()

#\, \mathrm{s}^{-1} \, \mathrm{cm}^{-2} \, \mathrm{Hz}^{-1} \, \mathrm{ster}^{-1})
