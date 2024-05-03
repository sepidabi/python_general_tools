import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
#import spectral as s
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec


# DECLERATIONS
indx_case = 1
fibril_indxs = ['06_142705','06_165307','06_165721','07_130930'] # ['06_165307', '06_165307', '06_165307', '06_165307']
bg_indxs =  ['06_142809','06_165404','06_165752','07_131005'] # ['06_165404', '06_165404', '06_165404', '06_165404']
index = '06_165307'
index_bg = '06_165404'
indxs = fibril_indxs
res = 0.0375 # CHROMIS pixel size in arcse
#tp = [2, 35, 75] # pixels of interest along the fibril
#gg = len(tp) # number of pixels of interest

# plotting info
plot_mode = 1 # 1 for each fibril profile
                          # 2 for (manual) dI histogram
                          # 3 for seaborn histograms
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }
markers = ['s','o','^'] # for the desired pixels
alp = 0.008 #transparency of the plotted lines
alp2 = 0.9


datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
savedir = datadir+'OUTPUT/'
resdir = '/scratch/sepid/stic_old/'
outdir = '/home/seki2695/OUTPUT/inv/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
fr = 28 #frame no.
ss = 0 # stokes param index
#xx = 100  # set to zero to choose the first pixel
#yy = 0     #stokes I values in the dat array

#calibration info
calib_cak = mf.readfits(datadir+'calib.3934.fits')
calib_ha = mf.readfits(datadir+'calib.6563.fits')

# fibril slab
fibdir = datadir+'fr'+str(fr)+'/'
fib_file = (file_search(fibdir,'crispex*3950*'+index+'.csav'))[0]
fib = restore(fibdir+fib_file)
f_slab = fib.loop_slab[:,fr,:]*1000  #intensity values in the desired frame
f_x = fib.x_coords
f_x_pts = fib.x_loop_pts
f_y = fib.y_coords
f_l = fib.loop_size
f_y_pts = fib.y_loop_pts


#fibril in Ha
fib_file_h = (file_search(fibdir,'crispex*6563*'+index+'.csav'))[0]
fib_h = restore(fibdir+fib_file_h)
f_slab_h = fib_h.loop_slab[:,fr,:]/calib_ha[0]  #intensity values in the desired frame

# bg slab
bg_file = file_search(fibdir,'crispex*3950*'+index_bg+'.csav')[0]
bg = restore(fibdir+bg_file)
b_slab = bg.loop_slab[:,fr,:]*1000 #intensity values in the desired frame
b_x = bg.x_coords
b_x_pts = bg.x_loop_pts
b_y_pts = bg.y_loop_pts
b_y = bg.y_coords
b_l = bg.loop_size

# fibril/background valid length
f_n_valid = np.min([f_l,b_l])
n = 20
tp = np.array([f_n_valid/n, f_n_valid/2, (n-1)*f_n_valid/n], dtype='uint')
gg = len(tp)


# bg in Ha
bg_file_h = (file_search(fibdir,'crispex*6563*'+index_bg+'.csav'))[0]
bg_h = restore(fibdir+bg_file_h)
b_slab_h = bg_h.loop_slab[:,fr,:]/calib_ha[0]  #intensity values in the desired frame

# To get the CaK spectral positions
file_ck = savedir+sorted(file_search(savedir, 'f_obs3950_'+index+'.fits'))[0]
cak_spect = mf.readfits(file_ck)[0,:,0]

# FIGURE
# setting the graphics window
plt.close('all')
f = plt.figure(figsize = [4,10])

left = 0.13
right = 0.87
wspace = 0
bottom = 0.07
top = 0.99
hspace = 0

gs = gridspec.GridSpec(1, 2) # grid scale
gs.update(left=left,
          right=right,
          wspace = wspace,
          bottom=bottom,
          top=top,
          hspace=hspace)

ax1 = host_subplot(gs[0,0], adjustable='box')
ax2 = host_subplot(gs[0,1], adjustable='box')

inrange = 5
plt_indx = np.linspace(5, 225, f_n_valid//inrange +1, dtype=int)
space = 1
ymin = 1e1*f_slab[0,plt_indx[-1]]/f_slab[-1,plt_indx[-1]]+ (f_n_valid - len(plt_indx))
ymax = 1e1*f_slab[0,plt_indx[0]]/f_slab[-1,plt_indx[0]]+ (f_n_valid)
ymin_b = 1e1*b_slab[20,plt_indx[-1]]/b_slab[-1,plt_indx[-1]]+ (f_n_valid - len(plt_indx))
ymax_b = 1e1*b_slab[20,plt_indx[0]]/b_slab[-1,plt_indx[0]]+ (f_n_valid)

for i in range(len(plt_indx)):

    x_end_points = np.array([cak_spect[20]-cak_spect[10], cak_spect[-1]-cak_spect[10]])
    y_end_points = 1e1*np.array([f_slab[20, plt_indx[i]]/f_slab[-1,plt_indx[i]] , b_slab[0, plt_indx[i]]/b_slab[-1,plt_indx[i]]])+ (f_n_valid - i*space)

    #ax1.plot(cak_spect[0:21]-cak_spect[10], 1e1*f_slab[0:21, plt_indx[i]]/f_slab[-1,plt_indx[i]] + (f_n_valid - i*space), alpha = 1.-plt_indx[i]*0.004, color = 'r')
    ax1.plot(cak_spect[0:21]-cak_spect[10], 1e1*f_slab[0:21, plt_indx[i]]/f_slab[-1,plt_indx[i]] + (f_n_valid - i*space), color = 'r')
    #print(1e1*b_slab[21, plt_indx[i]]/b_slab[-1,plt_indx[i]]+ (f_n_valid - i*space))
    #ax2.plot(cak_spect[0:21]-cak_spect[10], 1e1*b_slab[0:21, plt_indx[i]]/b_slab[-1,plt_indx[i]] + (f_n_valid - i*space), alpha = 1.-plt_indx[i]*0.004, color = 'Grey')
    ax2.plot(cak_spect[0:21]-cak_spect[10], 1e1*b_slab[0:21, plt_indx[i]]/b_slab[-1,plt_indx[i]] + (f_n_valid - i*space), color = 'Grey')
    #ax1.plot(x_end_points, y_end_points)
    
ax1.set_xlabel('$\Delta\lambda$ $\mathrm{[\AA]}$', fontsize = 8)
ax2.set_xlabel('$\Delta\lambda$ $\mathrm{[\AA]}$', fontsize = 8)
#ax1.set_xlim([1,19])
ax1.set_ylim([ymin,ymax+1])
ax1.set_yticks([ymin+1,np.mean([ymin,ymax]),ymax])
#ax1.set_xticks([10])
#ax1.set_xticklabels([0])
ax1.set_yticklabels([ur'$\u25b2$', ur'$\u25cf$', ur'$\u25a0$'])
ax1.set_ylabel(r'$I\mathrm{_f}$ / $I\mathrm{_{cont}}$', fontsize = 8)
ax1.tick_params(labelsize = 8, axis = 'both')
#ax2.set_xlim([1,19])
ax2.set_ylim([ymin,ymax+1])
ax2.set_yticks([ymin_b+1,np.mean([ymin_b,ymax_b]),ymax_b])
ax2.set_yticklabels([ur'$\u25b2$', ur'$\u25cf$', ur'$\u25a0$'])
ax2.yaxis.tick_right()
#ax2.set_xticks([10])
#ax2.set_xticklabels([0])
ax2.set_ylabel('$I\mathrm{_b}$ / $I\mathrm{_{cont}}$', fontsize = 8)
ax2.yaxis.set_label_position("right")
ax2.tick_params(labelsize = 8, axis = 'both')

plt.show()
filename = outdir+'cak_intensity_vary_nofade.pdf'
f.savefig(filename, resolution = 1000)
