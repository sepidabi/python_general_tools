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
#file_ca8 = file_search(datadir,'crispex*'+pref[1]+'*.fcube')
#file_cak = file_search(datadir,'crispex*'+pref[2]+'*.fcube')
#file_ha = file_search(datadir,'crispex*'+pref[3]+'*.fcube')
#cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])
#cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])
#cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])
#cube_ha = lp_read(datadir+file_ha[0],datadir+file_ha[1])
#c_map = cube_cak[fr,0,-1,:,:]
#cak_int_un = (lp_read_scan(datadir+'cak_int_un.fcube'))[fr,:,:]
#ha_int_un = (lp_read_scan(datadir+'h_int_un.fcube'))[fr,:,:]
w1_ca8, w2_ca8 = 7, 13
#ca8_int_un = np.mean(cube_ca8[fr, 0, w1_ca8:w2_ca8+1, :,:], axis = 0)

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


# Graphics
# figure sizes
dx, dy = 8., 5.
dx_grid = 5.
dx1, dx2 = dx*2./dx_grid, dx*3./dx_grid
plt.close('all')
f2, axs = plt.subplots(4,3,figsize=(dx2,dy))

plt.subplots_adjust(left = 0.1,
                    right = 0.97,
                    bottom = 0.09,
                    top = 0.95,
                    wspace = 0.27,
                    hspace = 0.20
)




#for fn in range(len(cat_indx)):
for fn in range(len(file_search(fibdir,'crispex*3950*2018*.csav'))/2):

    int_amp = 1.
    
    # extracting FIBRIL
    # Ca K
    fib_file = (file_search(fibdir,'crispex*3950*2018*.csav'))[2*fn]#[cat_indx[fn]]
    fib = restore(fibdir+fib_file)
    f_slab = fib.loop_slab[:,fr,:]*1000  # intensity values in the desired frame
    f_Itot = int_amp*np.mean(f_slab[cak_i:cak_n,:]/f_slab[-1,:], axis = 1)*100
    # Ha
    fib_file_h = (file_search(fibdir,'crispex*6563*2018*.csav'))[2*fn]#[cat_indx[fn]]
    fib_h = restore(fibdir+fib_file_h)
    f_slab_h = fib_h.loop_slab[:,fr,:]/calib_ha[0]  # intensity values in the desired frame
    f_Itot_h = int_amp*np.mean(f_slab_h[ha_i:ha_n,:]/f_slab_h[0,:], axis = 1)*100
    # Ca 8542
    fib_file_ca8 = (file_search(fibdir,'crispex*8542*2018*.csav'))[2*fn]#[cat_indx[fn]]
    fib_ca8 = restore(fibdir+fib_file_ca8)
    f_slab_ca8 = fib_ca8.loop_slab[0,:,fr,:]/calib_ca8[0]  # intensity values in the desired frame
    f_Itot_ca8 = int_amp*np.mean(f_slab_ca8[ca8_i:ca8_n,:]/f_slab_ca8[0,:], axis = 1)*100
    # coords
    f_x = fib.x_coords
    f_x_pts = fib.x_loop_pts
    f_y = fib.y_coords
    f_l = fib.loop_size
    f_y_pts = fib.y_loop_pts

    # extracting BG
    bg_file = file_search(fibdir,'crispex*3950*.csav')[2*fn+1]#[cat_indx[fn]+1]
    bg = restore(fibdir+bg_file)
    b_slab = bg.loop_slab[:,fr,:]*1000 #intensity values in the desired frame
    b_Itot = int_amp*np.mean(b_slab[cak_i:cak_n,:]/b_slab[-1,:], axis = 1)*100
    # Ha
    bg_file_h = (file_search(fibdir,'crispex*6563*.csav'))[2*fn+1]#[cat_indx[fn]+1]
    bg_h = restore(fibdir+bg_file_h)
    b_slab_h = bg_h.loop_slab[:,fr,:]/calib_ha[0]  #intensity values in the desired frame
    b_Itot_h = int_amp*np.mean(b_slab_h[ha_i:ha_n,:]/b_slab_h[0,:], axis = 1)*100
    # Ca 8542
    bg_file_ca8 = (file_search(fibdir,'crispex*8542*.csav'))[2*fn+1]#[cat_indx[fn]+1]
    bg_ca8 = restore(fibdir+bg_file_ca8)
    b_slab_ca8 = bg_ca8.loop_slab[0,:,fr,:]/calib_ca8[0]  #intensity values in the desired frame
    b_Itot_ca8 = int_amp*np.mean(b_slab_ca8[ca8_i:ca8_n,:]/b_slab_ca8[0,:], axis = 1)*100
    
    # coords
    b_x = bg.x_coords
    b_x_pts = bg.x_loop_pts
    b_y_pts = bg.y_loop_pts
    b_y = bg.y_coords
    b_l = bg.loop_size

    if (2*fn==cat_indx[0] or 2*fn==cat_indx[1] or 2*fn==cat_indx[2] or 2*fn==cat_indx[3]):

        cat_lab = ['f1','f2','f3','f4'] # category labels
        colors = ['coral', 'lightsalmon', 'lightgrey','darkgrey'] # category label colors
        lwidth = 2. # width of the lines
        alp = 0.3 #transparency of the plotted lines

        for ii in range (len(cat_indx)):
            if (2*fn==cat_indx[ii]):
                # intensity profiles
                ax11 = axs[ii,0]
                ax12 = axs[ii,1]
                ax13 = axs[ii,2]
                ax11.plot(cak_wav[cak_i:cak_n], f_Itot, color = 'red')
                ax11.plot(cak_wav[cak_i:cak_n], b_Itot, color = 'gray')
                ax12.plot(ha_wav[ha_i:ha_n], f_Itot_h, color = 'red')
                ax12.plot(ha_wav[ha_i:ha_n], b_Itot_h, color = 'gray')
                ax13.plot(ca8_wav[ca8_i:ca8_n], f_Itot_ca8, color = 'red')
                ax13.plot(ca8_wav[ca8_i:ca8_n], b_Itot_ca8, color = 'gray')

                # xaxis ticks
                xtick_pos_cak = [-1, -0.5, 0, 0.5, 1.]
                xtick_lab_cak = np.array(np.round(xtick_pos_cak), dtype = int)
                xtick_pos_ha = [-0.8, -0.4, 0, 0.4, 0.8]
                xtick_lab_ha = np.array(np.round(xtick_pos_ha), dtype = int)
                xtick_pos_ca8 = [-0.7, -0.3, 0, 0.3, 0.7]
                xtick_lab_ca8 = np.array(np.round(xtick_pos_ca8), dtype = int)
                cak_xlim = [-1.0,1.0]
                ha_xlim = [-0.8,0.8]
                ca8_xlim = [-0.7,0.7]
                ax11.set_xticks([])
                ax11.set_xlim(cak_xlim)
                ax13.set_ylabel(cat_lab[ii],fontdict = font)
                ax13.yaxis.set_label_position("right")
                ax12.set_xticks([])
                ax12.set_xlim(ha_xlim)
                ax13.set_xlim(ca8_xlim)
                ax13.set_xticks([])
                ax11.tick_params(axis ='both', labelsize =8)
                ax12.tick_params(axis ='both', labelsize =8)
                ax13.tick_params(axis ='both', labelsize =8)
 
                if ii==3:
                    ax11.set_xticks(xtick_pos_cak)
                    ax12.set_xticks(xtick_pos_ha)
                    ax13.set_xticks(xtick_pos_ca8)
                    #ax11.set_xlabel(r'$\lambda-\lambda_{0}$ $\mathrm{[\AA]}$', fontdict = font)
                    ax12.set_xlabel(r'$\lambda-\lambda_{0}$ $\mathrm{[\AA]}$', fontdict = font)
                    #ax13.set_xlabel(r'$\lambda-\lambda_{0}$ $\AA$', fontdict = font)
                    ax11.tick_params(axis ='both', labelsize =8)
                    ax12.tick_params(axis ='both', labelsize =8)
                    ax13.tick_params(axis ='both', labelsize =8)
                if ii==0:
                    ax11.set_title(r'Ca II K', fontdict = font)
                    ax12.set_title(r'H$\mathrm{\alpha}$', fontdict = font)
                    ax13.set_title(r'Ca II 8542 $\mathrm{\AA}$', fontdict = font)
                    
                # plot labels
                #x_lab_pos, y_lab_pos = 0, 3.65
                
                # lambda-integration range
                ax11.axvline(x = cak_wav[6], linestyle = ':', color = 'black', linewidth = 0.5)
                ax11.axvline(x = cak_wav[14], linestyle = ':', color = 'black', linewidth = 0.5)
                ax12.axvline(x = ha_wav[4], linestyle = ':', color = 'black', linewidth = 0.5)
                ax12.axvline(x = ha_wav[10], linestyle = ':', color = 'black', linewidth = 0.5)
                ax13.axvline(x = ca8_wav[7], linestyle = ':', color = 'black', linewidth = 0.5)
                ax13.axvline(x = ca8_wav[13], linestyle = ':', color = 'black', linewidth = 0.5)
                    
                # yaxis ticks
                cak_ylim = [8,18.5]
                ha_ylim = [16.5,82]
                ca8_ylim = [33,81]
                ax11.set_ylim(cak_ylim)
                ax12.set_ylim(ha_ylim)
                ax13.set_ylim(ca8_ylim)
                
f2.text(0.015,0.55,r'$\bar{I}$ / $I_\mathrm{cont}$ [%]', fontdict=font, rotation = 90)#f2.set_tight_layout(True)
plt.show()
filename = outdir+'category2.pdf'
plt.savefig(filename, quality = 100)
print 'file save to: '+ filename
