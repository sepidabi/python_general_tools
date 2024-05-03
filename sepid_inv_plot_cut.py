from __future__ import print_function
import sparsetools as sp
import spectral as s
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from Tkinter import *
import Tkinter
from scipy import ndimage
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from sepid import *

# plotting font, etc.
plot_mode = 3    # set to 1 for fibril prof
                             # 2 for RF
                             # 3 for along the fibril
                             # 4 for extra radiation plot
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }
markers = ['s','o','^'] # for the desired pixels

# settings
mode = 'roinew_fullstokes'
cyc = '1'
dir = '/home/seki2695/INV/stic/fov/'
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
savedir = datadir+'OUTPUT/'
outdir = '/home/seki2695/OUTPUT/inv/'
amp=3.5 # cak profile amplifier
fr = 28 # frame number of the data
# ROI coords
xmin, xmax = 1064, 1353
ymin, ymax = 487, 639
# CUT coords
x11, x12, x13, x14 = 80, 115, 189, 224
y11, y12, y13, y14 = 6, 42, 110, 147
rotate = 0
fmark = 0

# extract inversion results
i = sp.profile(dir+'roinew_observed_06_165307.nc')
i.dat[0,:,:,0:42,0]  = i.dat[0,:,:,0:42,0]*amp
o = sp.profile(dir+mode+'_synthetic_cycle1_06_165307.nc')
o.dat[0,:,:,0:42,0]  = o.dat[0,:,:,0:42,0]*amp
m = sp.model(dir+mode+'_atmosout_cycle1_06_165307.nc')

# fibril calibrated files
ffile_fe = file_search(savedir, 'f_obs6302_*.fits')
bfile_fe = file_search(savedir, 'b_obs6302_*.fits')

# specify the inversion cubes
temp = m.temp[0,:,:,:].squeeze()*1e-3
vlos = m.vlos[0,:,:,:].squeeze()*1e-5
vmt = m.vturb[0,:,:,:].squeeze()*1e-5

# cube dimensions
ny, nx, ndep = temp.shape
dep = m.ltau[0,0,0,:]
# Specifying the depth in which to check the results
dep_n = 28

# perpendicular cut -> coords has to be manually chosen
x1, y1, x2, y2 = 200, 46, 170, 132
slop = (y2-y1) / (x2-x1) # cut slope
# rotation angle
theta = 90 - math.degrees(np.arctan(abs(slop)))

# overplotting the slope
#if rotate==0:
#    ax.plot([x1,x2], [y1,y2], color = 'black', linestyle = '--')

# rotating cubes
temp_rot = ndimage.rotate(temp, theta, reshape=True)
vlos_rot = ndimage.rotate(vlos, theta, reshape=True)
vmt_rot = ndimage.rotate(vmt, theta, reshape=True)
ny_rot, nx_rot, ndep = temp_rot.shape

# manually chosen coords on the rotated map
x1_rot, x2_rot = 110, 225
y1_rot, y2_rot = 60, 170

# define projection matrix 1
unrot = np.ones((ny,nx,ndep), dtype = float)

# each fibril process
for fn in range(len(ffile_fe)):
    fibril = mf.readfits(savedir + ffile_fe[fn])
    dark = mf.readfits(savedir + bfile_fe[fn])
            
    # fibril slab
    fibdir = datadir+'fr'+str(fr)+'/'
    fib_file = (file_search(fibdir,'crispex*3950*.csav'))[2*fn]
    fib = restore(fibdir+fib_file)
    f_slab = fib.loop_slab[:,fr,:]*1000  #intensity values in the desired frame
    f_x = fib.x_coords
    f_x_pts = fib.x_loop_pts
    f_y = fib.y_coords
    f_l = fib.loop_size
    f_y_pts = fib.y_loop_pts

    # marking the fibril coords
    if fmark==1:
        for i in range (len(f_x)):
            # check if the coords are inside the ROI
            fx, fy = int(f_x[i]-xmin), int(f_y[i]-ymin)
            if fx > 0 and fx < nx -1 and fy > 0 and fy < ny-1:
                # setting the coord to ZERO in the projection matrix
                unrot[fy, fx,:] = 0.

    # bg slab
    bg_file = file_search(fibdir,'crispex*3950*.csav')[2*fn+1]
    bg = restore(fibdir+bg_file)
    b_slab = bg.loop_slab[:,fr,:]*1000 #intensity values in the desired frame
    b_x = bg.x_coords
    b_x_pts = bg.x_loop_pts
    b_y_pts = bg.y_loop_pts
    b_y = bg.y_coords
    b_l = bg.loop_size

    
    # marking the bg coords
    if fmark ==0:
        for i in range (len(b_x)):
            bx, by = int(b_x[i]-xmin), int(b_y[i]-ymin)
            if bx > 0 and bx < nx -1 and by > 0 and by < ny-1:
                # setting the coord to ZERO in the projection matrix
                unrot[by, bx,:] = 0

    # fibril and bg overplots
    lwidth = 0.5
    alp = 0.5
    if rotate==0:
        
        plt.imshow(temp[:,:,dep_n], vmax = 7)
        
        plt.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'grey', alpha = alp, linewidth = lwidth)
        plt.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'red', alpha = alp, linewidth = lwidth)

# rotated projection matrix
rot = ndimage.rotate(unrot, theta, reshape=True)

# marked temperature map
finale = np.multiply(rot, temp_rot)[:,:,dep_n]

# clipping the depth between log = [0,-7]
dep_init = 5 # logtau = -7
dep_fin = 53 # logtau = 0
dep_cropped = dep[dep_init: dep_fin+1]

##########################
#                         Head Cut
##########################
'''
!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! BE WARE !!!!!! The following indexes have been detected manually
!!!!!!!!!!!!!!!!!!!!!!!! according to fibrils location on the rotated ROI map
'''
tp_f_head = np.array([73, 113, 124, 127, 135, 142, 163]) - y1_rot
tp_b_head = np.array([68, 121, 122, 131, 139, 145, 158]) - y1_rot

tp_f_tail = np.array([84, 112, 123, 165]) - y1_rot
tp_b_tail = np.array([75, 118, 120, 158]) - y1_rot

for n in range(2):
    # head cut
    if n==0:
        # y values
        tp_f = tp_f_head
        tp_b = tp_b_head
        # x value
        xcut = x2_rot
        markers = ['o', 's', '*', 'P', '^', 'D', 'p']
        filename = 'inv_perpendicular_headcut.pdf'
    if n==1:
        # y values
        tp_f = tp_f_tail
        tp_b = tp_b_tail
        # x value
        xcut = x1_rot
        markers = ['o', 's', '*', 'P']
        filename = 'inv_perpendicular_tailcut.pdf'

    # cut profile maps
    temperature = np.transpose(temp_rot[y1_rot:y2_rot, xcut, dep_init:dep_fin+1])
    vlos = np.transpose(vlos_rot[y1_rot:y2_rot, xcut, dep_init:dep_fin+1])
    vmt = np.transpose(vmt_rot[y1_rot:y2_rot, xcut, dep_init:dep_fin+1])
        
    # graphics
    plt.close('all')

    f = plt.figure(figsize = [5.5,6])
    gs = gridspec.GridSpec(3, 1) # grid scale of the 1st col.
    aspect = 1 # panels aspect ratio (y/x)
    gs.update(left=0.075,
              right=0.88,
              #wspace=0.05,
              bottom=0.08,
              top=0.99,
              hspace = 0.0)
        
    ax1 = plt.subplot(gs[0,0],adjustable = 'box',aspect = 'equal')
    ax2 = plt.subplot(gs[1,0],adjustable = 'box',aspect = 'equal')
    ax3 = plt.subplot(gs[2,0],adjustable = 'box',aspect = 'equal')
                
    # RF depth grids
    depth_grid_sp = 7
    ytick_pos = np.linspace(0,6,7,dtype='uint')*depth_grid_sp
    ytick_lab = np.round(dep_cropped[ytick_pos])
        
    # temperature
    tmax = 9
    tmin = 4
    panel1 = ax1.imshow(temperature, cmap = 'afmhot', vmin = tmin, vmax = tmax,aspect = aspect)

    ax1.set_ylabel(r'log($\mathrm{\tau_{500}}$)', fontdict=font)
    ax1.set_xlim(0, temperature.shape[1]-1)
    ax1.set_xticklabels([])
    ax1.set_yticks(ytick_pos)
    ax1.set_yticklabels(ytick_lab.astype(int), fontdict=font)
    for pp in range(len(tp_f)):

        lwidth = 1.
        alph = 0.65
        
        # fibril indicator
        ax1.axvline(tp_f[pp], color = 'red', linestyle = '--', linewidth = lwidth, alpha = alph)
        ax2.axvline(tp_f[pp], linestyle = '--', color = 'red', linewidth = lwidth, alpha = alph)
        ax3.axvline(tp_f[pp], linestyle = '--', color = 'red', linewidth = lwidth, alpha = alph)
        # fibril marker
        #ax1.plot(tp_f[pp],3, marker = markers[pp], markersize = 5, color = 'red')

        # bg indocator
        ax1.axvline(tp_b[pp], color = 'grey', linestyle = '--', linewidth = lwidth, alpha = alph)
        ax2.axvline(tp_b[pp], linestyle = '--', color = 'grey', linewidth = lwidth, alpha = alph)
        ax3.axvline(tp_b[pp], linestyle = '--', color = 'grey', linewidth = lwidth, alpha = alph)
        # fibril marker
        #ax1.plot(tp_b[pp],3, marker = markers[pp], markersize = 5, color = 'blue')
        
    # temp colorbar
    axins = inset_axes(ax1, #axis
                       width="2%",  # width = 10% of parent_bbox width
                       height="90%",  # height : 50%
                       loc='center left',
                       bbox_to_anchor=(1.01, 0., 1, 1),
                       bbox_transform=ax1.transAxes,
                       borderpad=0,
    )
    ct_ticks = np.linspace(tmin,tmax+1,tmax-tmin+1,dtype='uint')
    ct = plt.colorbar(panel1, cax=axins, orientation="vertical",
                      ticks = ct_ticks,
    )
    ct.set_label(r'T (kK)', fontdict = font)
    ct.ax.set_yticklabels(['4', '5', '6', '7', '8', '9', '10', '$>$11'], fontdict=font)
        
    # vlos
    norm_v = np.max(np.abs(vlos))
    panel2 = ax2.imshow(vlos, cmap = 'bwr', vmin = -norm_v, vmax = norm_v,aspect = aspect)
    ax2.set_ylabel(r'log($\mathrm{\tau_{500}}$)', fontdict=font)
    ax2.set_xlim(0, vlos.shape[1]-1)
    ax2.set_xticklabels([])
    ax2.set_yticks(ytick_pos)
    ax2.set_yticklabels(ytick_lab.astype(int), fontdict=font)
        
    # temp colorbar
    axins = inset_axes(ax2, #axis
                       width="2%",  # width = 10% of parent_bbox width
                       height="90%",  # height : 50%
                       loc='center left',
                       bbox_to_anchor=(1.01, 0., 1, 1),
                       bbox_transform=ax2.transAxes,
                       borderpad=0,
    )
    cv = plt.colorbar(panel2, cax=axins, orientation="vertical")
    cv.set_label(r'v$_{\rm LOS}$ [km/s]', fontdict=font)
        
    # vmt
    norm_vmt = np.max(np.abs(vmt))
    panel3 = ax3.imshow(vmt, cmap = 'bone', vmin = 0, vmax = norm_vmt,aspect = aspect)
    ax3.set_xlabel(r'pixel along the cut', fontdict=font)
    ax3.set_xlim(0, temperature.shape[1]-1)
    ax3.set_ylabel(r'log($\mathrm{\tau_{500}}$)', fontdict=font)
    ax3.set_yticks(ytick_pos)
    ax3.set_yticklabels(ytick_lab.astype(int), fontdict=font)
        
    # vmt colorbar
    axins = inset_axes(ax3, #axis
                       width="2%",  # width = 10% of parent_bbox width
                       height="90%",  # height : 50%
                       loc='center left',
                       bbox_to_anchor=(1.01, 0., 1, 1),
                       bbox_transform=ax3.transAxes,
                       borderpad=0,
    )
    cvmt = plt.colorbar(panel3, cax=axins, orientation="vertical")
    cvmt.set_label(r'v$_{\rm turb}$ [km/s]', fontdict = font)
        
    plt.show()

    outname = outdir+filename
    f.savefig(outname, quality = 100)
    print ('\nFile saved to: ' + outname)
    #break

temp_cube = np.flip(temp_rot[y1_rot:y2_rot, x1_rot:x2_rot, :], axis = 2)
np.save(outdir+'cool_cube2.npy', temp_cube)

# define projection matrix 2
backrot  = np.ones((ny_rot,nx_rot, ndep), dtype = float)
backrot[y1_rot,x1_rot,:] = 0
backrot[y1_rot,x2_rot,:] = 0
backrot[y2_rot,x1_rot,:] = 0
backrot[y2_rot,x2_rot,:] = 0
finale = np.multiply(backrot,temp_rot)

# cropping coords
x1_crop, x2_crop = 45, 334
y1_crop, y2_crop = 86, 238

# rotating back and cropping edges
check_backrot = ndimage.rotate(finale, -theta, reshape=True)[y1_crop:y2_crop, x1_crop:x2_crop,:]

# CUT coords
x11, x12, x13, x14 = 80, 115, 189, 224
y11, y12, y13, y14 = 6, 42, 110, 147

# side cuts
plt.plot([x11,x13], [y13,y14], color = 'white')
plt.plot([x12,x14], [y11,y12], color = 'white')

# perpen. cuts
plt.plot([x12,x11], [y11,y13], color = 'white', linestyle = '--')
plt.plot([x14,x13], [y12,y14], color = 'white', linestyle = '--')

# TO VISUALIZE THE CUBE in home laptop

'''
# in terminal:
ETS_TOOLKIT=wx
ipython --gui=wx

# in ipyton
import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.interactive(True)

# packages
from mayavi import mlab
import numpy as np

# load the cube
cub = np.load('/Users/seki2695/Desktop/Fibril/cool_cube2.npy')
s = cub[:,:,7:47] # between 0<log(tau)<-6
y,x,z = s.shape
s = np.flip(s, axis=0)

# visualize
b1 = np.percentile(s, 20)
b2 = np.percentile(s, 80)
mlab.pipeline.volume(mlab.pipeline.scalar_field(s), vmin=b1, vmax=b2)
mlab.axes(ranges=[0,y-1, 0, x-1, 0, -6], nb_labels = 7, zlabel = 'log(tau)', xlabel = '[pxls]', ylabel = '[pxls]')
mlab.colorbar(title = 'T [kK]', orientation = 'vertical')
mlab.show()

'''
