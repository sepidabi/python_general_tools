# sepid_plot_category.py
'''
The script is meant to demonstrate a subFOV
with sample of categorised bright fibrils overplotted
'''

# import modules
import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
#import spectral as s

# DECLERATIONS
cat_indx = [202, 196, 190, 182]#[182,190,196,202] # fibril count index in crispex
xmin, xmax = 1064, 1353 #999+50+20, 1423-20-20-10
ymin, ymax = 487, 639 #457+25+20-30,744-10-70-25 # cropping the ROI

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

# directories paths
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
invdir = '/scratch/sepid/stic_old/fov/'
savedir = datadir+'OUTPUT/'
resdir = '/scratch/sepid/stic/'
outdir = '/home/seki2695/OUTPUT/inv/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
fr = 28 #frame no.
fibdir = datadir+'fr'+str(fr)+'/'
ss = 0 # stokes param index
res = 0.0375 # CHROMIS pixel size in arcse
dep = 28 # RFmax depth

i = sp.profile(invdir+'roinew_observed_06_165307.nc')
i.dat[0,:,:,0:42,0]  = i.dat[0,:,:,0:42,0]#*amp
o = sp.profile(invdir+'roinew_smoothed_synthetic_cycle2_06_165307.nc')
o.dat[0,:,:,0:42,0]  = o.dat[0,:,:,0:42,0]#*amp
m = sp.model(invdir+'roinew_smoothed_atmosout_cycle2_06_165307.nc')
temp = m.temp[0,:,:,dep]

# maps from the data
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
#file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
file_ca8 = file_search(datadir,'crispex*'+pref[1]+'*.fcube')
file_cak = file_search(datadir,'crispex*'+pref[2]+'*.fcube')
file_ha = file_search(datadir,'crispex*'+pref[3]+'*.fcube')
#cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])
cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])
cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])
cube_ha = lp_read(datadir+file_ha[0],datadir+file_ha[1])
c_map = cube_cak[fr,0,-1,:,:]
cak_int_un = (lp_read_scan(datadir+'cak_int_un.fcube'))[fr,:,:]
ha_int_un = (lp_read_scan(datadir+'h_int_un.fcube'))[fr,:,:]
w1_ca8, w2_ca8 = 7, 13
ca8_int_un = np.mean(cube_ca8[fr, 0, w1_ca8:w2_ca8+1, :,:], axis = 0)

# fibril calibrated files
#ffile_fe = file_search(savedir, 'f_obs6302_*.fits')
#bfile_fe = file_search(savedir, 'b_obs6302_*.fits')

# calibration info and wave settings
calib_cak = mf.readfits(datadir+'calib.3934.fits')
calib_ha = mf.readfits(datadir+'calib.6563.fits')
calib_ca8 = mf.readfits(datadir+'calib.8542.fits')
ha_spec_pos = mf.readfits(datadir+'ha_wave.fits')
ca8_spec_pos = mf.readfits(datadir+'ca8_wave.fits')
ha_c = 6562.8
ha_wav = ha_spec_pos# * 1.e3 # -> mAA
ca8_wav = ca8_spec_pos# * 1.e3 # -> mAA
cak_spec = mf.readfits(datadir+'cak_wave.fits')
cak_wav = (cak_spec - cak_spec[10]) *1.e10 #1.e13 # ->AA -> mAA
# wave length crop
cak_i, cak_n = 0,21
ha_i, ha_n = 0, 15
ca8_i, ca8_n = 0, -1


# ROI maps
cak_map_cropped = cak_int_un[ymin:ymax,xmin:xmax].squeeze()*10.e7
ha_map_cropped = ha_int_un[ymin:ymax,xmin:xmax].squeeze()/10.e3
ca8_map_cropped = ca8_int_un[ymin:ymax,xmin:xmax].squeeze()/10.e3

# Graphics
# figure sizes
dx, dy = 8., 5.
dx_grid = 5.
dx1, dx2 = dx*2./dx_grid, dx*3./dx_grid
plt.close('all')
f1 = plt.figure(figsize = (dx1,dy))
ax1 = f1.add_subplot(3,1,1)
ax2 = f1.add_subplot(3,1,2)
ax3 = f1.add_subplot(3,1,3)

plt.subplots_adjust(left = 0.13,
                    right = 0.99999,
                    bottom = 0.08,
                    top = 0.97,
                    wspace = 0.0,
                    hspace = 0.06
)

# maps
ax1.imshow(unsharp(cak_map_cropped), cmap = 'gray', aspect = 'equal', clim = (0.,5.5))
ax2.imshow(unsharp(ha_map_cropped), cmap = 'gray', aspect = 'equal')
ax3.imshow(unsharp(unsharp(ca8_map_cropped)), cmap = 'gray', aspect = 'equal', vmin = 0.0068, vmax = 0.0188)


# maps axes
ax1.set_xlim([0, xmax-xmin])
ax1.set_ylim([0, ymax-ymin])
ax2.set_ylabel('y [arcsec]', fontdict = font)
ax3.set_ylabel('y [arcsec]', fontdict = font)
ax1.set_xlabel('', fontdict = font)
ytick_pos = np.arange(0,(np.round((ymax-ymin)/40)+1)*40,40)
ytick_lab = ytick_pos*res
xtick_pos = np.arange(0,(np.round((xmax-xmin)/40)+1)*40,40)
xtick_lab = xtick_pos*res
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
ax1.tick_params(which='minor', length=2)
ax1.set_yticks(ytick_pos)
ax1.set_yticklabels(ytick_lab, fontdict=font)
ax1.set_xticks(xtick_pos)
ax1.set_xticklabels([])
ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
ax2.tick_params(which='minor', length=2)
ax2.set_xlim(0, xmax-xmin)
ax2.set_ylim(0, ymax-ymin)
ax1.set_ylabel('y [arcsec]', fontdict = font)
ax2.set_yticks(ytick_pos)
ax2.set_yticklabels(ytick_lab, fontdict=font)
ax2.set_xlabel('', fontdict = font)
ax2.set_xticks(xtick_pos)
ax2.set_xticklabels([])
ax3.xaxis.set_minor_locator(AutoMinorLocator(5))
ax3.yaxis.set_minor_locator(AutoMinorLocator(5))
ax3.tick_params(which='minor', length=2)
ax3.set_xlim(0, xmax-xmin)
ax3.set_ylim(0, ymax-ymin)
ax3.set_yticks(ytick_pos)
ax3.set_yticklabels(ytick_lab, fontdict=font)
ax3.set_xlabel('x [arcsec]', fontdict = font)
ax3.set_xticks(xtick_pos)
ax3.set_xticklabels([])
ax3.set_xticklabels(xtick_lab, fontdict=font)





# PQRS BOX
# cut coords
x11, x12, x13, x14 = 80, 115, 189, 224
y11, y12, y13, y14 = 6, 42, 110, 147

# side cuts
ax1.plot([x11,x13], [y13,y14], color = 'white', linewidth = 1)
ax1.plot([x12,x14], [y11,y12], color = 'white', linewidth = 1)
#ax2.plot([x11,x13], [y13,y14], color = 'white', linewidth = 1)
#ax2.plot([x12,x14], [y11,y12], color = 'white', linewidth = 1)
#ax3.plot([x11,x13], [y13,y14], color = 'white', linewidth = 1)
#ax3.plot([x12,x14], [y11,y12], color = 'white', linewidth = 1)

# perpen. cuts
ax1.plot([x12,x11], [y11,y13], color = 'white', linestyle = '--', linewidth = 1)
ax1.plot([x14,x13], [y12,y14], color = 'white', linestyle = '--', linewidth = 1)
#ax2.plot([x12,x11], [y11,y13], color = 'white', linestyle = '--', linewidth = 1)
#ax2.plot([x14,x13], [y12,y14], color = 'white', linestyle = '--', linewidth = 1)
#ax3.plot([x12,x11], [y11,y13], color = 'white', linestyle = '--', linewidth = 1)
#ax3.plot([x14,x13], [y12,y14], color = 'white', linestyle = '--', linewidth = 1)

# point labels
ax1.text(x14+2, y12-2, 'P', color = 'white', fontsize = 7)
ax1.text(x13+2, y14-2, 'Q', color = 'white', fontsize = 7)
ax1.text(x12-7, y11-6, 'R', color = 'white', fontsize = 7)
ax1.text(x11-7, y13, 'S', color = 'white', fontsize = 7)



# Vertical axis label
#f1.text(0.001,0.46,r'I$_{mean}$ [10$^{-5}$ erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$]', fontdict=font, rotation = 90)
# Horizontal Axis label
#f1.text(0.5,0.005,r'$\mathrm{\Delta \lambda}$ [$\mathrm{\AA}$]', fontdict=font)
ax1.text(1,2, 'Ca II K', color = 'white', fontdict = font, alpha = 0.7)
ax2.text(1,2,  'H$\mathrm{\\alpha}$', color = 'white', fontdict = font, alpha = 0.7)
ax3.text(1,2,  'Ca II 8542 $\mathrm{\AA}$', color = 'white', fontdict = font, alpha = 0.7)


#for fn in range(len(cat_indx)):
for fn in range(len(file_search(fibdir,'crispex*3950*.csav'))/2):

    int_amp = 1.
    
    # extracting FIBRIL
    # Ca K
    fib_file = (file_search(fibdir,'crispex*3950*.csav'))[2*fn]#[cat_indx[fn]]
    fib = restore(fibdir+fib_file)
    f_slab = fib.loop_slab[:,fr,:]*1000  # intensity values in the desired frame
    f_Itot = int_amp*np.mean(f_slab[cak_i:cak_n,:]/f_slab[-1,:], axis = 1)*100
    # Ha
    fib_file_h = (file_search(fibdir,'crispex*6563*.csav'))[2*fn]#[cat_indx[fn]]
    fib_h = restore(fibdir+fib_file_h)
    f_slab_h = fib_h.loop_slab[:,fr,:]/calib_ha[0]  # intensity values in the desired frame
    f_Itot_h = int_amp*np.mean(f_slab_h[ha_i:ha_n,:]/f_slab_h[0,:], axis = 1)*100
    # coords
    f_x = fib.x_coords
    f_x_pts = fib.x_loop_pts
    f_y = fib.y_coords
    f_l = fib.loop_size
    f_y_pts = fib.y_loop_pts

    # extracting BG
    bg_file = file_search(fibdir,'crispex*.csav')[2*fn+1]#[cat_indx[fn]+1]
    bg = restore(fibdir+bg_file)
    b_slab = bg.loop_slab[:,fr,:]*1000 #intensity values in the desired frame
    b_Itot = int_amp*np.mean(b_slab[cak_i:cak_n,:]/b_slab[-1,:], axis = 1)*100
    # Ha
    bg_file_h = (file_search(fibdir,'crispex*6563*.csav'))[2*fn+1]#[cat_indx[fn]+1]
    bg_h = restore(fibdir+bg_file_h)
    b_slab_h = bg_h.loop_slab[:,fr,:]/calib_ha[0]  #intensity values in the desired frame
    b_Itot_h = int_amp*np.mean(b_slab_h[ha_i:ha_n,:]/b_slab_h[0,:], axis = 1)*100
    
    # coords
    b_x = bg.x_coords
    b_x_pts = bg.x_loop_pts
    b_y_pts = bg.y_loop_pts
    b_y = bg.y_coords
    b_l = bg.loop_size

    # plot labels
    x_lab_pos, y_lab_pos = 0, 3.65

    if (2*fn==cat_indx[0] or 2*fn==cat_indx[1] or 2*fn==cat_indx[2] or 2*fn==cat_indx[3]):

        cat_lab = ['f1','f2','f3','f4'] # category labels
        colors = ['coral', 'lightsalmon', 'lightgrey','darkgrey'] # category label colors
        lwidth = 2. # width of the lines
        alp = 0.6 #transparency of the plotted lines

        ax1.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'orangered', linewidth = lwidth, alpha = alp, linestyle = '--')
        ax1.plot(b_x_pts-xmin,b_y_pts-ymin, color = 'black', linewidth = lwidth, alpha = alp, linestyle = '--')
        ax2.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'orangered', linewidth = lwidth, alpha = alp, linestyle = '--')
        ax3.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'orangered', linewidth = lwidth, alpha = alp, linestyle = '--')

        # intensity profiles
        alp2 = 0.8
        if (2*fn == cat_indx[0]):
            # sample labels
            ax1.text(f_x_pts[-1]-xmin-18, f_y_pts[-1]-ymin, cat_lab[0], alpha = alp2, color = 'orangered')
        

            
        if (2*fn == cat_indx[1]):
            # sample labels
            ax1.text(f_x_pts[-1]-xmin-18, f_y_pts[-1]-ymin, cat_lab[1], alpha = alp2, color = 'orangered')
            
        if (2*fn == cat_indx[2]):
            # sample labels
            ax1.text(f_x_pts[-1]-xmin-18, f_y_pts[-1]-ymin, cat_lab[2], alpha = alp2, color = 'orangered')
            
        if (2*fn == cat_indx[3]):
            # sample labels
            ax1.text(f_x_pts[-1]-xmin-18, f_y_pts[-1]-ymin, cat_lab[3], alpha = alp2, color = 'orangered')
                
    else:
        lwidth = 0.85 # width of the lines
        alp = 0.4 #transparency of the plotted lines

        ax1.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'yellow', linewidth = lwidth, alpha = alp, linestyle = '-')
        ax2.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'yellow', linewidth = lwidth, alpha = alp, linestyle = '-')
        ax3.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'yellow', linewidth = lwidth, alpha = alp, linestyle = '-')

#f1.set_tight_layout(True)
plt.show()
filename = outdir+'category1_modified.pdf'
plt.savefig(filename, quality = 100)
print 'file save to: '+ filename
