from __future__ import print_function
import sparsetools as sp
import spectral as s
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


dir='/nadir-scratch/jaime/2016.09.19/inv_jorrit/'
dir='/home/sdani/nadir-scratch/sdani/chromis/20160919/loop/run/'
dir='/home/sdani/nadir-scratch/jaime/loop/run/'
dir='stic/example_me/'
i = sp.profile(dir+'b_observed_06_165307.nc')
o = sp.profile(dir+'b_synthetic_cycle1_06_165307.nc')
m = sp.model(dir+'b_atmosout_cycle1_06_165307.nc')
#o = sp.profile(dir+'synthetic_cycle2_reg350p.nc')
#m = sp.model(dir+'atmosout_cycle2_reg350p.nc')
print(i.dat.shape)
print(o.dat.shape)
print(m.ltau.shape)
nt, nx, ny, nw, ns = i.dat.shape
#  chi2=0.0
#  npchi=0
#for ww in range(25,62,1):
#for ww in range(112,147,1):    
#    if o.dat[0,0,0,ww,0] > 0.2:
#        chi2+=(i.dat[0,:,:,ww,0]-o.dat[0,:,:,ww,0])**2
#          npchi+=1
#  for ww in range(178,217,1):
#      if o.dat[0,0,0,ww,0] > 0.2:
#          chi2+=(i.dat[0,:,:,ww,0]-o.dat[0,:,:,ww,0])**2
#          npchi+=1

#chi2/=(npchi*ns)
#print(npchi)

fig = plt.figure()
ax = plt.subplot2grid((3,3), (0,0), colspan=2, rowspan=2)
ax2 = plt.subplot2grid((3,3), (2,0), colspan=2, rowspan=1)
ax3 = plt.subplot2grid((3,3), (0,2), colspan=1, rowspan=1)
ax4 = plt.subplot2grid((3,3), (1,2), colspan=1, rowspan=1)
ax5 = plt.subplot2grid((3,3), (2,2), colspan=1, rowspan=1)
logt=[20,23,35]
#-4.5
print(m.ltau[0,0,0,logt[0]])
#ax.imshow(m.temp[0,:,:,logt[0]]/1.e3,vmax=10.,aspect=0.1)
#ax.imshow(m.B[0,:,:,logt[0]]/1.e3,vmax=10.)
ax.imshow(m.vlos[0,:,:,logt[0]]*1.e-5,vmax=10.,vmin=-10,aspect=0.1)
#ax.imshow(m.pgas[0,:,:,0])
#ax.imshow(m.inc[0,:,:,logt[0]])
#ax.imshow(chi2,vmax=50.)



def onclick(event):
    xx, yy = int(event.xdata), int(event.ydata)
    print('xdata=%d, ydata=%d' %
          (xx, yy))
    ax2.cla()
    ax2.plot(i.dat[0,yy,xx,:,0],'.', color='black')
    ax2.plot(o.dat[0,yy,xx,:,0], color='orangered')
    ax2.set_ylabel('Intensity')
    ax2.set_xlabel(r'Wavelength [$\mathrm{\AA}$]')
    #ax2.set_ylim(0.05, 3)

    ax3.cla()   
    ax3.plot(m.ltau[0,yy,xx,:],m.temp[0,yy,xx,:]/1.e3, 'k-')
    ax3.set_ylabel('T [kK]')


    ax4.cla()   
    ax4.plot(m.ltau[0,yy,xx,:],m.vturb[0,yy,xx,:]*1.e-5, 'k-')
    ax4.set_ylabel('v$_{turb}$ [km/s]')
    #ax4.plot(m.ltau[0,yy,xx,:],np.log(m.pgas[0,yy,xx,:]), 'k-')
    #ax4.set_ylabel('Pg [dyn/cm2]')

    ax5.cla()   
    ax5.plot(m.ltau[0,yy,xx,:],m.vlos[0,yy,xx,:]*1.e-5, 'k-')
    ax5.set_ylabel('v$_{los}$ [km/s]')
    ax5.set_xlabel('ltau')
    
    fig.canvas.draw()
    


    
cid = fig.canvas.mpl_connect('button_press_event', onclick)


plt.show()

