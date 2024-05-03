# This script is meant to plot
# the inversion results of the FOV

import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import seaborn as sns
import scipy.ndimage as spnd
import matplotlib
from matplotlib.gridspec import GridSpec


# fonts
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }

#DECLERATIONS
index = '06_165307'
index_bg = '06_165404'
xmin, xmax, ymin, ymax = 245, 1678, 10, 1130 # FOV coords

# specify the directory and files
resdir = '/home/seki2695/INV/stic/fov/'
dir = resdir + 'results_nocmap/'
outdir = '/home/seki2695/OUTPUT/inv/'
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
savedir = datadir+'OUTPUT/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']

# data extract
res = 0.0375 # CHROMIS pixel size in arcsec
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
fr = 28 #frame no.
file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
file_ca8 = file_search(datadir,'crispex*'+pref[1]+'*.fcube')
file_cak = file_search(datadir,'crispex*'+pref[2]+'*.fcube')
file_ha = file_search(datadir,'crispex*'+pref[3]+'*.fcube')
cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])[fr]
cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])[fr]
cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])[fr]
cube_ha = lp_read(datadir+file_ha[0],datadir+file_ha[1])[fr]
c_map = cube_cak[0,-1,:,:]
cak_int_un = (lp_read_scan(datadir+'cak_int_un.fcube'))[fr,:,:]

# read the fov data from its fits file
header = fits.getheader(resdir + 'fov_inv.fits')
print(header)
fov = fits.getdata(resdir + 'fov_inv.fits')
fov_cmap = np.transpose(fits.getdata(resdir + 'fov_backup_no_cmap.fits'))

# dimensions
yy, xx, zz = fov.shape[0], fov.shape[1], fov.shape[2]

# tau grid
tau = fov[0,0,:,-1]

# fibril calibrated files
ffile_fe = file_search(savedir, 'f_obs6302_*.fits')
bfile_fe = file_search(savedir, 'b_obs6302_*.fits')

# maps
fov_temp = fov[:,:,:,0]
fov_vlos = fov[:,:,:,1]
fov_vlos_cmap = fov_cmap[:,:,:,1]
fov_vturb = fov[:,:,:,2]
fov_Bln = fov[:,:,:,3]
fov_Bln_mod = fov[:,:,:,8]

#'''
#########
# FOV figures
#########
depths = [51, 23]
dep = depths[1]
plt.close('all')
f = plt.figure(figsize = [7.2,8.75-1.75])
gs = gridspec.GridSpec(3,2) # grid scale of the 1st col.
aspect =1 # dunno!
gs.update(left=0.025,
          right=0.93,
          wspace=0.025,
          bottom=0.051,
          top=0.97,
          hspace = 0.03,
)

# axis info
sc_fact = 266.7 # axis scaling factor
ytick_pos = np.arange(0,(np.round(yy/sc_fact)+1)*sc_fact,sc_fact)
ytick_lab = np.round(ytick_pos*res).astype(int)
xtick_pos = np.arange(0,(np.round(xx/sc_fact)+1)*sc_fact,sc_fact)
xtick_lab = np.round(xtick_pos*res).astype(int)


for ii in range(len(depths)):
    dd = depths[ii]
    # axis settings
    ax_temp = plt.subplot(gs[0,ii],adjustable = 'box')
    ax_temp.set_xticks([])
    ax_temp.set_xlim(0,xx)
    ax_temp.set_ylim(0,yy)
    ax_temp.set_yticks(ytick_pos)
    ax_temp.set_ylabel(r'y [arcsec]', fontdict=font)
    ax_temp.set_yticklabels(ytick_lab, fontdict = font)
    ax_temp.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax_temp.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax_temp.tick_params(which='minor', length=2)
    ax_temp.set_xticks(xtick_pos)
    ax_temp.set_xticklabels('')
    
    ax_vlos = plt.subplot(gs[1,ii],adjustable = 'box')#,aspect = 'equal')
    ax_vlos.set_xticks([])
    ax_vlos.set_xlim(0,xx)
    ax_vlos.set_ylim(0,yy)
    ax_vlos.set_yticks(ytick_pos)
    ax_vlos.set_yticklabels(ytick_lab, fontdict = font)
    ax_vlos.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax_vlos.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax_vlos.tick_params(which='minor', length=2)
    ax_vlos.set_xticks(xtick_pos)
    ax_vlos.set_xticklabels('')
    ax_vlos.set_ylabel(r'y [arcsec]', fontdict=font)

    ax_B = plt.subplot(gs[2,ii],adjustable = 'box')#,aspect = 'equal')
    ax_B.set_xticks(xtick_pos)
    ax_B.set_xticklabels([])
    ax_B.set_xlim(0,xx)
    ax_B.set_ylim(0,yy)
    ax_B.set_yticks(ytick_pos)
    ax_B.set_yticklabels(ytick_lab, fontdict = font)
    ax_B.set_ylabel(r'y [arcsec]', fontdict=font)
    ax_B.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax_B.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax_B.tick_params(which='minor', length=2)
    ax_B.set_xticks(xtick_pos)
    ax_B.set_xlabel(r'x [arcsec]', fontdict=font)
    ax_B.set_xticklabels(xtick_lab, fontdict = font)

    #ax_vturb = plt.subplot(gs[3,ii],adjustable = 'box')#,aspect = 'equal')
    #ax_vturb.set_xticks(xtick_pos)
   # ax_vturb.set_xticklabels(xtick_lab, fontdict = font)
   # ax_vturb.set_xlim(0,xx)
   # ax_vturb.set_ylim(0,yy)
   # ax_vturb.set_yticks(ytick_pos)
   # ax_vturb.set_yticklabels(ytick_lab, fontdict = font)
   # ax_vturb.set_xlabel(r'x [arcsec]', fontdict=font)
   # ax_vturb.xaxis.set_minor_locator(AutoMinorLocator(10))
   # ax_vturb.yaxis.set_minor_locator(AutoMinorLocator(10))
   # ax_vturb.tick_params(which='minor', length=2)
   # ax_vturb.set_xticks(xtick_pos)
   # ax_vturb.set_xlabel(r'x [arcsec]', fontdict=font)
    #ax_vturb.set_ylabel(r'y [arcsec]', fontdict=font)

        
    # temperature maps
    temp_min, temp_max = np.min(fov_temp[:, :, dd]), np.max(fov_temp[:, :, dd])
    if ii==1:
        temp_min, temp_max = 4.7, 7.5
    else:
        temp_min, temp_max = 5.5, 7.
    panel_temp = ax_temp.imshow(unsharp(fov_temp[:,:,dd], alpha = 0.5, sigma = 2.),
                                cmap = 'gist_heat', origin = 'lower', vmin = temp_min, vmax = temp_max)
    ax_temp.set_title(r'log($\tau_{\rm {500}}$ ) = '+str(np.round(tau[dd], decimals = 2)), fontdict = font)

    # vlos map
    if ii ==0:
        vlos_norm = 5
    else:
        vlos_norm = 8
    panel_vlos = ax_vlos.imshow(fov_vlos[:, :, dd],
                                cmap = 'bwr', vmin = -vlos_norm, vmax = vlos_norm, origin = 'lower')

    # B map
    B_norm = np.max(abs(fov_Bln[:, :, dd]))
    
    panel_B = ax_B.imshow(fov_Bln[:, :, dd], cmap = 'RdGy_r', vmin = -1.5, vmax = 1.5, origin = 'lower')

    # vturb map
    #vturb_norm = np.max(abs(fov_vturb[:,:,dd]))
    #panel_vturb = ax_vturb.imshow(np.transpose(fov_vturb[:,:,dd]), cmap = 'bone', vmin = 0, vmax = vturb_norm, origin = 'lower')

    # temp colorbar
    axin_temp = inset_axes(ax_temp,
                           width="3%",  # width = 10% of parent_bbox width
                           height="90%",  # height : 50%
                           loc='center left',
                           bbox_to_anchor=(1.02, 0., 1, 1),
                           bbox_transform=ax_temp.transAxes,
                           borderpad=0,
    )
    #ct_ticks =np.arange(np.ceil(temp_min), np.floor(temp_max), 0.5)
    #ct_tick_labels = 
    ct = plt.colorbar(panel_temp, cax=axin_temp, orientation="vertical")#, ticks = ct_ticks)
    ct.ax.tick_params(labelsize=8) 
    
    # vlos colorbar
    axin_vlos = inset_axes(ax_vlos, #axis
                           width="3%",  # width = 10% of parent_bbox width
                           height="90%",  # height : 50%
                           loc='center left',
                           bbox_to_anchor=(1.02, 0., 1, 1),
                           bbox_transform=ax_vlos.transAxes,
                           borderpad=0,
    )
    cv = plt.colorbar(panel_vlos, cax=axin_vlos, orientation="vertical")
    cv.ax.tick_params(labelsize=8) 


    # B colorbar axis
    axin_B = inset_axes(ax_B, #axis
                   width="3%",  # width = 10% of parent_bbox width
                        height="90%",  # height : 50%
                        loc='center left',
                        bbox_to_anchor=(1.02, 0., 1, 1),
                        bbox_transform=ax_B.transAxes,
                        borderpad=0,
    )
    cB = plt.colorbar(panel_B, cax=axin_B, orientation="vertical")
    cB.ax.tick_params(labelsize=8) 

    # vmt colorbar axis
    #axin_vturb = inset_axes(ax_vturb, #axis
      #                      width="3%",  # width = 10% of parent_bbox width
        #                    height="90%",  # height : 50%
          #                  loc='center left',
            #                bbox_to_anchor=(1.02, 0., 1, 1),
              #              bbox_transform=ax_vturb.transAxes,
                #            borderpad=0,
    #)
    #cvmt = plt.colorbar(panel_vturb, cax=axin_vturb, orientation="vertical")
    #cvmt.ax.tick_params(labelsize=8) 

    if ii == 1:
        #ax_vturb.set_yticklabels([])
        ax_vlos.set_yticklabels([])
        ax_temp.set_yticklabels([])
        ax_B.set_yticklabels([])
        ax_B.set_ylabel('')
        ax_temp.set_ylabel('')
        ax_vlos.set_ylabel('')
        #ax_vturb.set_ylabel('')
        ct.set_label(r'$T$ [kK]', fontdict = font)
        cv.set_label(r'$v_{\rm LOS}$ [km s${\rm^{-1}}$]', fontdict=font)
        cB.set_label(r'$B_{\rm long}$ [kG]', fontdict = font)
        #cvmt.set_label(r'v$_{\rm turb}$ [km s${\rm^{-1}}$]', fontdict = font)

        # fibrils overplot
        lwidth = 0.4 # width of the lines
        alp = 0.85 #transparency of the plotted lines

        # extracting ALL fibrilS
        n_pxls = 2351 # the number was set by hand 
        fb_temp = np.zeros(0)
        bg_temp = np.zeros(0)
        fb_vlos = np.zeros(0)
        bg_vlos = np.zeros(0)
        fb_vturb = np.zeros(0)
        bg_vturb= np.zeros(0)

        for fn in range(len(ffile_fe)):
            #print('still working...')
            # Ca K
            fibdir = datadir+'fr'+str(fr)+'/'
            fib_file = (file_search(fibdir,'crispex*3950*.csav'))[2*fn]#[cat_indx[fn]]
            bg_file = (file_search(fibdir,'crispex*3950*.csav'))[2*fn+1]#[cat_indx[fn]]
            fb = restore(fibdir+fib_file)
            bg = restore(fibdir+bg_file)
            # fib coords
            f_x = fb.x_coords
            f_x_pts = fb.x_loop_pts
            f_y = fb.y_coords
            f_l = fb.loop_size
            f_y_pts = fb.y_loop_pts
            # bg coords
            b_x = bg.x_coords
            b_x_pts = bg.x_loop_pts
            b_y = bg.y_coords
            b_l = bg.loop_size
            b_y_pts = bg.y_loop_pts


            ax_temp.plot(f_x_pts-xmin,f_y_pts-ymin,linestyle = ':', color = 'white', linewidth = lwidth, alpha = alp)
            ax_vlos.plot(f_x_pts-xmin,f_y_pts-ymin,linestyle = ':', color = 'black', linewidth = lwidth, alpha = alp)
            ax_B.plot(f_x_pts-xmin,f_y_pts-ymin,linestyle = ':', color = 'black', linewidth = lwidth, alpha = alp)
            #ax_vturb.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'white', linewidth = lwidth, alpha = alp)

            for j in range(np.min([len(f_x), len(b_x)])):
                if (np.int(f_y[j]) < ymax and np.int(f_x[j]) <xmax and
                    np.int(b_y[j]) < ymax and np.int(b_x[j]) <xmax):
                    fb_temp = np.append(fb_temp, fov_temp[np.int(f_y[j])-ymin, np.int(f_x[j])-xmin, :])
                    fb_vlos = np.append(fb_vlos, fov_vlos_cmap[np.int(f_y[j])-ymin, np.int(f_x[j])-xmin, :]-0.5)
                    fb_vturb = np.append(fb_vturb, fov_vturb[np.int(f_y[j])-ymin, np.int(f_x[j])-xmin, :])
                    bg_temp = np.append(bg_temp, fov_temp[np.int(b_y[j])-ymin, np.int(b_x[j])-xmin, :])
                    bg_vlos = np.append(bg_vlos, fov_vlos_cmap[np.int(b_y[j])-ymin, np.int(b_x[j])-xmin, :]-0.5)
                    bg_vturb = np.append(bg_vturb, fov_vturb[np.int(b_y[j])-ymin, np.int(b_x[j])-xmin, :])

# making the cubes of the fibrillar and background pixels along the depths
pxls_ndep = fb_temp.size
ndep = len(tau) # number of tau grids
pxls = pxls_ndep/ndep
# reshaping them...
fb_temp = fb_temp.reshape((pxls,ndep))
fb_vlos = fb_vlos.reshape((pxls,ndep))
fb_vturb = fb_vturb.reshape((pxls,ndep))
bg_temp = bg_temp.reshape((pxls,ndep))
bg_vlos = bg_vlos.reshape((pxls,ndep))
bg_vturb = bg_vturb.reshape((pxls,ndep))

plt.show()
filename = outdir + 'test_fov_fullstokes_maps.pdf'
f.savefig(filename, dpi = 1000)
print 'file saved to: '+ filename
plt.show()


##############
#         PROFILES        #
##############

# Visual Statistics
plt.close('all')
alp = 0.5

f = plt.figure(figsize=(7.5,3.))

# fonts
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 7.,
        }

# figure margins
left = 0.08
bottom = 0.11
right = 0.99
top = 0.98
wspace = 0.86
hspace = 0.33

ax21 = plt.subplot2grid((3,6), (0,0), colspan=2, rowspan=3)
ax22 = plt.subplot2grid((3,6), (0,2), colspan=2, rowspan=3)
ax23 = plt.subplot2grid((3,6), (0,4), colspan=2, rowspan=3)
#ax11 = plt.subplot2grid((3,6), (2,0), colspan=2, rowspan=1)
#ax12 = plt.subplot2grid((3,6), (2,2), colspan=2, rowspan=1)
#ax13 = plt.subplot2grid((3,6), (2,4), colspan=2, rowspan=1)

# means and sigmas
fb_temp_m = np.mean(fb_temp, axis = 0)
fb_vlos_m = np.mean(fb_vlos, axis = 0)
fb_vturb_m = np.mean(fb_vturb, axis = 0)
fb_temp_s = np.std(fb_temp, axis = 0)
fb_vlos_s = np.std(fb_vlos, axis = 0)
fb_vturb_s = np.std(fb_vturb, axis = 0)
bg_temp_m = np.mean(bg_temp, axis = 0)
bg_vlos_m = np.mean(bg_vlos, axis = 0)
bg_vturb_m = np.mean(bg_vturb, axis = 0)
bg_temp_s = np.std(bg_temp, axis = 0)
bg_vlos_s = np.std(bg_vlos, axis = 0)
bg_vturb_s = np.std(bg_vturb, axis = 0)

alpha, lw= 0.002, 0.5
lw_s, alpha_s = 0.5, 0.2
lw_m, alpha_m = 0.8, 0.9

# temp
#for jj in range(pxls):
 #   ax21.plot(bg_temp[jj,:], tau, color = 'grey', alpha = alpha, linewidth = lw)
  #  ax21.plot(fb_temp[jj,:], tau, color = 'red', alpha = alpha, linewidth = lw)
ax21.fill_betweenx(tau, bg_temp_m - bg_temp_s, bg_temp_m + bg_temp_s, alpha = alpha_s, color = 'grey')
ax21.fill_betweenx(tau, fb_temp_m - fb_temp_s, fb_temp_m + fb_temp_s, alpha = alpha_s, color = 'orangered')
ax21.plot(bg_temp_m, tau, color = 'black', alpha = alpha_m, linewidth = lw_m, label = r'f')
ax21.plot(fb_temp_m, tau, color = 'orangered', alpha = alpha_m, linewidth = lw_m, label = r'b')
ax21.axhline(tau[dep], linestyle = ':', color = 'black', alpha = 0.5)
ax21.set_ylim([-2,-5.5])
ax21.set_xlim([4.2,7.8])
ax21.set_xlabel(r'$T$ [kK]', fontdict = font)
ax21.set_ylabel(r'log($\tau_\mathrm{{500}}$)', fontdict = font)
ax21.tick_params(labelsize = 7)


# vlos
#for jj in range(pxls):
 #   ax22.plot(bg_vlos[jj,:], tau, color = 'grey', alpha = alpha, linewidth = lw)
  #  ax22.plot(fb_vlos[jj,:], tau, color = 'red', alpha = alpha, linewidth = lw)
ax22.fill_betweenx(tau, bg_vlos_m - bg_vlos_s, bg_vlos_m + bg_vlos_s, alpha = alpha_s, color = 'grey')
ax22.fill_betweenx(tau, fb_vlos_m - fb_vlos_s, fb_vlos_m + fb_vlos_s, alpha = alpha_s, color = 'orangered')
ax22.plot(bg_vlos_m, tau, color = 'black', alpha = alpha_m, linewidth = lw_m)
ax22.plot(fb_vlos_m, tau, color = 'orangered', alpha = alpha_m, linewidth = lw_m)
ax22.axhline(tau[dep], linestyle = ':', color = 'black', alpha = 0.5)
ax22.set_ylim([-2,-5.5])
ax22.set_xlim([-3.5,3.5])
ax22.set_xlabel(r'$v_\mathrm{{LOS}}$ [km s${\rm^{-1}}$]', fontdict = font)
ax22.tick_params(labelsize = 7)

# vturb
#for jj in range(pxls):
 #   ax23.plot(bg_vturb[jj,:], tau, color = 'grey', alpha = alpha, linewidth = lw)
  #  ax23.plot(fb_vturb[jj,:], tau, color = 'red', alpha = alpha, linewidth = lw)
ax23.fill_betweenx(tau, bg_vturb_m - bg_vturb_s, bg_vturb_m + bg_vturb_s, alpha = alpha_s, color = 'grey')
ax23.fill_betweenx(tau, fb_vturb_m - fb_vturb_s, fb_vturb_m + fb_vturb_s, alpha = alpha_s, color = 'orangered')
ax23.plot(bg_vturb_m, tau, color = 'black', alpha = alpha_m, linewidth = lw_m)
ax23.plot(fb_vturb_m, tau, color = 'orangered', alpha = alpha_m, linewidth = lw_m)
ax23.axhline(tau[dep], linestyle = ':', color = 'black', alpha = 0.5)
ax23.set_ylim([-2,-5.5])
ax23.set_xlim([0,6])
ax23.set_xlabel(r'$v_\mathrm{{turb}}$ [km s${\rm^{-1}}$]', fontdict = font)
ax23.tick_params(labelsize = 7)
#f.legend(loc = 'upper center', ncol = 2, fontsize = 8, markerscale = 0.1, frameon = False, columnspacing = 1., handlelength = 1., handletextpad = 0.3)

matplotlib.pyplot.subplots_adjust(left=0.08
                                  , bottom=0.12
                                  , right=0.99
                                  , top=0.98
                                  , wspace=0.86
                                  , hspace=0.33
)
#f.tight_layout()
filename = outdir + 'test_fov_fullstokes_stats.pdf'
f.savefig(filename, dpi = 1000)
print 'file saved to: '+ filename
plt.show()

# joint plot figure settings

# figure size in inches
fsize_x = 3 
fsize_y = 3

# figure margins
left = 0.20
bottom = 0.17
right = 0.99
top = 0.99
wspace = 0.11
hspace = 0.11

fontsize = 9.

################
#     TEMPERATURE
################

plt.close('all')
fig = plt.figure(figsize = (fsize_x,fsize_y))

x = fb_temp[:,dep]
y = bg_temp[:,dep]

xmin, xmax = temp_min, 7

gs = GridSpec(5,5)

ax_joint = fig.add_subplot(gs[1:5,0:4])
ax_marg_x = fig.add_subplot(gs[0,0:4])
ax_marg_y = fig.add_subplot(gs[1:5,4])

# Turn off ticks and labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
plt.setp(ax_marg_y.get_xticklabels(), visible=False)
plt.setp(ax_marg_x.get_yticklabels(), visible=False)
ax_marg_x.get_yaxis().set_visible(False)
ax_marg_y.get_xaxis().set_visible(False)


# horizontal hist
ax_marg_x.spines['top'].set_visible(False)
ax_marg_x.spines['right'].set_visible(False)
ax_marg_x.spines['left'].set_visible(False)
# vertical hist
ax_marg_y.spines['top'].set_visible(False)
ax_marg_y.spines['right'].set_visible(False)
ax_marg_y.spines['bottom'].set_visible(False)
# joint plot
ax_joint.spines['top'].set_visible(False)
ax_joint.spines['right'].set_visible(False)

# Set labels on joint
ax_joint.set_xlabel(r'$T_\mathrm{{f}}$ [kK]', fontdict = font)
ax_joint.set_ylabel(r'$T_\mathrm{{b}}$ [kK]', fontdict = font)
sns.kdeplot(x,y, cmap = 'Reds', ax = ax_joint, shade = True)
sns.kdeplot(x, ax = ax_marg_x, shade=True, color = 'orangered')
sns.kdeplot(y, ax = ax_marg_y, shade=True, color = 'orangered', vertical = True)

# Set ax limits on marginals
ax_joint.set_xlim([xmin,xmax])
ax_joint.set_ylim([xmin,xmax])
ax_marg_y.set_ylim([xmin,xmax])
ax_marg_x.set_xlim([xmin,xmax])
ax_joint.tick_params(labelsize = fontsize)
ax_marg_x.tick_params(labelsize = fontsize)
ax_marg_y.tick_params(labelsize = fontsize)

pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
print(np.round(stat.pearsonr(x, y)[0],decimals = 5))
pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
ax_joint.annotate('pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')
ax_joint.text(6.5,6.65, 'y = x', fontsize = fontsize-2, rotation = 45, color = 'darkred', va = 'bottom', ha = 'left')

# extras overplotted lines
xx = np.linspace(xmin, xmax, 10)
ax_joint.plot(xx, xx, linestyle = ':', color = 'darkred', linewidth = 0.75)
#ax_joint.plot(xx, stat.pearsonr(x, y)[0]*xx, linestyle = '-', color = 'black', linewidth = 0.75)


plt.tight_layout()
plt.subplots_adjust(left = left,
                    bottom = bottom,
                    right = right,
                    top = top,
                    wspace = wspace,
                    hspace = hspace
)
plt.show()
fig.savefig(outdir + 'T_joint_logt' + str(np.round(tau[dep], decimals = 1)) + '.pdf', dpi = 1000)


################
#     VLOS
################

plt.close('all')
fig = plt.figure(figsize = (fsize_x,fsize_y))

x = fb_vlos[:,dep]
y = bg_vlos[:,dep]

xmin, xmax = -4, 6

gs = GridSpec(5,5)

ax_joint = fig.add_subplot(gs[1:5,0:4])
ax_marg_x = fig.add_subplot(gs[0,0:4])
ax_marg_y = fig.add_subplot(gs[1:5,4])

# Turn off ticks and labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
plt.setp(ax_marg_y.get_xticklabels(), visible=False)
plt.setp(ax_marg_x.get_yticklabels(), visible=False)
ax_marg_x.get_yaxis().set_visible(False)
ax_marg_y.get_xaxis().set_visible(False)


# horizontal hist
ax_marg_x.spines['top'].set_visible(False)
ax_marg_x.spines['right'].set_visible(False)
ax_marg_x.spines['left'].set_visible(False)
# vertical hist
ax_marg_y.spines['top'].set_visible(False)
ax_marg_y.spines['right'].set_visible(False)
ax_marg_y.spines['bottom'].set_visible(False)
# joint plot
ax_joint.spines['top'].set_visible(False)
ax_joint.spines['right'].set_visible(False)

# Set labels on joint
ax_joint.set_xlabel(r'$v_\mathrm{{LOS,f}}$ [km s${\rm^{-1}}$]', fontdict = font)
ax_joint.set_ylabel(r'$v_\mathrm{{LOS,b}}$ [km s${\rm^{-1}}$]', fontdict = font)
sns.kdeplot(x,y, cmap = 'Reds', ax = ax_joint, shade = True)
sns.kdeplot(x, ax = ax_marg_x, shade=True, color = 'orangered')
sns.kdeplot(y, ax = ax_marg_y, shade=True, color = 'orangered', vertical = True)

# Set ax limits on marginals
ax_joint.set_xlim([xmin,xmax])
ax_joint.set_ylim([xmin,xmax])
ax_marg_y.set_ylim([xmin,xmax])
ax_marg_x.set_xlim([xmin,xmax])
ax_joint.tick_params(labelsize = fontsize)
ax_marg_x.tick_params(labelsize = fontsize)
ax_marg_y.tick_params(labelsize = fontsize)

pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
print(np.round(stat.pearsonr(x, y)[0],decimals = 5))
pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
ax_joint.annotate('pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')
ax_joint.text(4.2,4.7, 'y = x', fontsize = fontsize-2, rotation = 45, color = 'darkred', va = 'bottom', ha = 'left')

# extras overplotted lines
xx = np.linspace(xmin, xmax, 10)
ax_joint.plot(xx, xx, linestyle = ':', color = 'darkred', linewidth = 0.75)
#ax_joint.plot(xx, stat.pearsonr(x, y)[0]*xx, linestyle = '-', color = 'black', linewidth = 0.75)


plt.tight_layout()
plt.subplots_adjust(left = left,
                    bottom = bottom,
                    right = right,
                    top = top,
                    wspace = wspace,
                    hspace = hspace
)
plt.show()
fig.savefig(outdir + 'vlos_joint_logt' + str(np.round(tau[dep], decimals = 1)) + '.pdf', dpi = 1000)


################
#     Vturb
################

plt.close('all')
fig = plt.figure(figsize = (fsize_x,fsize_y))

x = fb_vturb[:,dep]
y = bg_vturb[:,dep]

xmin, xmax = -1, 9

gs = GridSpec(5,5)

ax_joint = fig.add_subplot(gs[1:5,0:4])
ax_marg_x = fig.add_subplot(gs[0,0:4])
ax_marg_y = fig.add_subplot(gs[1:5,4])

# Turn off ticks and labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
plt.setp(ax_marg_y.get_xticklabels(), visible=False)
plt.setp(ax_marg_x.get_yticklabels(), visible=False)
ax_marg_x.get_yaxis().set_visible(False)
ax_marg_y.get_xaxis().set_visible(False)


# horizontal hist
ax_marg_x.spines['top'].set_visible(False)
ax_marg_x.spines['right'].set_visible(False)
ax_marg_x.spines['left'].set_visible(False)
# vertical hist
ax_marg_y.spines['top'].set_visible(False)
ax_marg_y.spines['right'].set_visible(False)
ax_marg_y.spines['bottom'].set_visible(False)
# joint plot
ax_joint.spines['top'].set_visible(False)
ax_joint.spines['right'].set_visible(False)

# Set labels on joint
ax_joint.set_xlabel(r'$v_\mathrm{{turb,f}}$ [km s${\rm^{-1}}$]', fontdict = font)
ax_joint.set_ylabel(r'$v_\mathrm{{turb,b}}$ [km s${\rm^{-1}}$]', fontdict = font)
sns.kdeplot(x,y, cmap = 'Reds', ax = ax_joint, shade = True)
sns.kdeplot(x, ax = ax_marg_x, shade=True, color = 'orangered')
sns.kdeplot(y, ax = ax_marg_y, shade=True, color = 'orangered', vertical = True)

# Set ax limits on marginals
ax_joint.set_xlim([xmin,xmax])
ax_joint.set_ylim([xmin,xmax])
ax_marg_y.set_ylim([xmin,xmax])
ax_marg_x.set_xlim([xmin,xmax])
ax_joint.tick_params(labelsize = fontsize)
ax_marg_x.tick_params(labelsize = fontsize)
ax_marg_y.tick_params(labelsize = fontsize)

pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
print(np.round(stat.pearsonr(x, y)[0],decimals = 5))
pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
ax_joint.annotate('pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')
ax_joint.text(6.9,7.5, 'y = x', fontsize = fontsize-2, rotation = 45, color = 'darkred', va = 'bottom', ha = 'left')

# extras overplotted lines
xx = np.linspace(xmin, xmax, 10)
ax_joint.plot(xx, xx, linestyle = ':', color = 'darkred', linewidth = 0.75)
#ax_joint.plot(xx, stat.pearsonr(x, y)[0]*xx, linestyle = '-', color = 'black', linewidth = 0.75)


plt.tight_layout()
plt.subplots_adjust(left = left,
                    bottom = bottom,
                    right = right,
                    top = top,
                    wspace = wspace,
                    hspace = hspace
)
plt.show()
fig.savefig(outdir + 'vturb_joint_logt' + str(np.round(tau[dep], decimals = 1)) + '.pdf', dpi = 1000)





##########################
# TEMPERATURE DIFFERENCE HISTOGRAM
##########################

plt.close('all')

fontsize = 8.

# fonts
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': fontsize,
        }

temp_dif = (fb_temp-bg_temp)*1.e3
temp_dif_m = np.mean(temp_dif, axis = 0)
temp_dif_s = np.std(temp_dif, axis = 0)
fig = plt.figure(figsize = [4,4])

for ii in range(len(tau)):
    plt.scatter(temp_dif[:,ii], np.zeros(temp_dif.shape[0],dtype = float)+tau[ii]
                , marker = 's'
                , color = 'darkred'
                , s = 50
                , alpha = 0.003
                , edgecolors = 'none'
    )

# data range plot
lw = 1.
plt.plot(temp_dif_m,tau,color = 'darkred', linewidth = lw)
plt.plot(temp_dif_s+temp_dif_m,tau,color = 'darkred', linestyle = '--', linewidth = lw)
plt.plot(-temp_dif_s+temp_dif_m,tau,color = 'darkred', linestyle = '--', linewidth = lw)

plt.ylim(-1,-5.5)
plt.xlim(-600,1250)
plt.xlabel(r'T [K]', fontdict = font)
plt.ylabel(r'log($\tau_{\mathrm{500}}$)', fontdict = font)
plt.tick_params(axis = 'both', labelsize = fontsize)

left = 0.16  # the left side of the subplots of the figure
right = 0.96   # the right side of the subplots of the figure
bottom = 0.11  # the bottom of the subplots of the figure
top = 0.96     # the top of the subplots of the figure
wspace = 0.2  # the amount of width reserved for space between subplots,
              # expressed as a fraction of the average axis width
hspace = 0.2  # the amount of height reserved for space between subplots,


plt.subplots_adjust(left = left,
                    bottom = bottom,
                    right = right,
                    top = top,
                    wspace = wspace,
                    hspace = hspace
)
filename = outdir + 'delta_T.pdf'
plt.savefig(filename, dpi = 1000)









'''
##############
#     HISTOGRAMS      #
##############
dep = 21
alp = 0.35
# temp
bins = np.linspace(4.5,8, 50)
ax11.hist(fb_temp[:,dep], bins = bins, alpha = alp, color = 'orangered', histtype='stepfilled', label = 'fibril')#, x='X1', **kwargs)
ax11.hist(bg_temp[:,dep], bins, alpha = alp, color = 'gray', histtype='stepfilled', label = 'background')#, x='X2', **kwargs)
ax11.set_xlim([4.5,8])
ax11.set_ylabel(r'# of locations', fontdict=font)
ax11.tick_params(labelsize = 8)
#ax11.legend()

# vlos
bins = np.linspace(-7,7, 50)

ax12.hist(fb_vlos[:,dep], bins = bins, alpha = alp, color = 'orangered', histtype='stepfilled')#, x='X1', **kwargs)
ax12.hist(bg_vlos[:,dep], bins = bins, alpha = alp, color = 'gray', histtype='stepfilled')#, x='X2', **kwargs)
ax12.set_xlabel(r'v$_{\rm LOS}$ [km s${\rm^{-1}}$]', fontdict=font)
ax12.set_xlim([-7,7])
#ax12.set_ylabel(r'# of locations', fontdict=font)
ax12.tick_params(labelsize = 8)

# vturb
bins = np.linspace(0,8, 50)
ax13.hist(fb_vturb[:,dep], bins = bins, alpha = alp, color = 'orangered', histtype='stepfilled')#, x='X1', **kwargs)
ax13.hist(bg_vturb[:,dep], bins = bins, alpha = alp, color = 'gray', histtype='stepfilled')#, x='X2', **kwargs)
ax13.set_xlabel(r'v$_{\rm turb}$ [km s${\rm^{-1}}$]', fontdict=font)
ax13.set_xlim([0,8])
#ax13.set_ylabel(r'# of locations', fontdict=font)
ax13.tick_params(labelsize = 8)
'''

