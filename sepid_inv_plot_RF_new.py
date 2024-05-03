import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
#import spectral as s

# DECLERATIONS
###########
indx_case = 1
fibril_indxs = ['06_142705','06_165307','06_165721','07_130930']
bg_indxs = ['06_142809','06_165404','06_165752','07_131005']
index = fibril_indxs[indx_case]
index_bg = bg_indxs[indx_case]
tp = [2, 35, 75] # pixels of interest along the fibril
gg = len(tp) # number of pixels of interest

# plotting font, etc.
plot_mode = 2    # set to 1 for fibril prof
                             # 2 for RF
                             # 3 for along the fibril
                             # 4 for extra radiation plot
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }
markers = ['s','o','^'] # for the desired pixels
alp = 0.01 #transparency of the plotted lines

# inversion info
if plot_mode ==2 or plot_mode==3:
    mode = 'RF' # set to 'RF' for the response function (grid)
else:
    mode = ''
cyc = '1'
# Same settings as XXVI -> for actual inversion
nodes_temp = -6.88333, -5.96667, -5.23333, -4.31667, -3.4, -2.48333, -1.56667, -0.833333, 0.083333
nodes_vlos =  -6.8, -4.46666667, -2.13333333, 0.2
nodes_vturb = -6.0, -4, -1.9, -0.2
nodes_blong = 2
nodes_bhor = 1
nodes_azi = 1
#(Temp, Vlos, vturb, B, inc, azi, pgas_boundary)
regularization_type = 1,2,3,1,1,0,1
regularize = 1.
regularization_weights = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0

# settings board!
reg = 'n$_{T}$ = '+" ".join(str('%.1f'%(x)) for x in nodes_temp)+' \nn$_{v_{LOS}}$ = '+" ".join(str('%.1f'%(x)) for x in nodes_vlos)+' \nn$_{v_{turb}}$ = '+" ".join(str('%.1f'%(x)) for x in nodes_vturb)+'  \nreg$_{type}$ = '+" ".join(str(int(x)) for x in regularization_type)+'        reg = '+str('%.1f'%(regularize))+'       reg$_{weight}$ = '+" ".join(str('%.1f'%(x)) for x in regularization_weights)

# data and output dir
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
outdir = '/home/seki2695/OUTPUT/inv/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
fn = 28 #frame no.
ss = 0 # stokes param index
#xx = 100  # set to zero to choose the first pixel
#yy = 0     #stokes I values in the dat array

# extracting inversion results
# fibril
f_i = sp.profile('/home/seki2695/INV/stic/fb_back/f'+mode+'_observed_'+index+'.nc')
f_o = sp.profile('/home/seki2695/INV/stic/fb_back/f'+mode+'_synthetic_cycle'+cyc+'_'+index+'.nc')
f_m = sp.model('/home/seki2695/INV/stic/fb_back/f'+mode+'_atmosout_cycle'+cyc+'_'+index+'.nc')
f_mRF = sp.model('/home/seki2695/INV/stic/fb_back/f'+mode+'_atmosout_cycle'+cyc+'_'+index+'.nc')
f_rf = sp.profile('/home/seki2695/INV/stic/fb_back/fRF_synthetic_cycle2_'+index+'.nc')
# backgraound
b_i = sp.profile('/home/seki2695/INV/stic/bg_back/b'+mode+'_observed_'+index+'.nc')
b_o = sp.profile('/home/seki2695/INV/stic/bg_back/b'+mode+'_synthetic_cycle'+cyc+'_'+index+'.nc')
b_m = sp.model('/home/seki2695/INV/stic/bg_back/b'+mode+'_atmosout_cycle'+cyc+'_'+index+'.nc')
b_rf = sp.profile('/home/seki2695/INV/stic/bg_back/bRF_synthetic_cycle2_'+index+'.nc')
b_mRF = sp.model('/home/seki2695/INV/stic/bg_back/b'+mode+'_atmosout_cycle'+cyc+'_'+index+'.nc')


#Fibril extract
##========
amp = 1. #cak amplifier
#inversion
f_obs = f_i.dat[0,0,:,:,ss]
f_obs[:,0:42]  = f_i.dat[0,0,:,0:42,ss]*amp
f_syn = f_o.dat[0,0,:,:,ss]
f_syn[:,0:42]  = f_o.dat[0,0,:,0:42,ss]*amp
f_wav = f_i.wav
f_n = f_i.dat.shape[2]
f_dep = f_m.ltau[0,0,0,:]
f_depRF = f_mRF.ltau[0,0,0,:]
# slab
fibdir = datadir+'fr'+str(fn)+'/'
fib_file = file_search(fibdir,'crispex*3950*'+index+'.csav')
fib = restore(fibdir+fib_file[0])
f_slab = fib.loop_slab  #intensity values in the desired frame
f_x = fib.x_coords
f_x_pts = fib.x_loop_pts
f_y = fib.y_coords
f_l = fib.loop_size
f_y_pts = fib.y_loop_pts

#Background extract
#===========
# inversion
b_obs = b_i.dat[0,0,:,:,ss]
b_obs[:,0:42] = b_i.dat[0,0,:,0:42,ss]*amp
b_syn = b_o.dat[0,0,:,:,ss]
b_syn[:,0:42] = b_o.dat[0,0,:,0:42,ss]*amp
b_wav = b_i.wav
b_n =  b_i.dat.shape[2]
b_dep = b_m.ltau[0,0,0,:]
# slab
bg_file = file_search(fibdir,'crispex*3950*'+index_bg+'*.csav')
bg = restore(fibdir+bg_file[0])
b_slab = bg.loop_slab #intensity values in the desired frame
b_x = bg.x_coords
b_x_pts = bg.x_loop_pts
b_y_pts = bg.y_loop_pts
b_y = bg.y_coords
b_l = bg.loop_size

# fibril/background valid length
f_n_valid = np.min([f_n,b_n])
n = 20
tp = np.array([f_n_valid/n, f_n_valid/2, (n-1)*f_n_valid/n], dtype='uint')
gg = len(tp)

#containers:
f_temp_m = np.zeros(45)
f_temp_med = np.zeros(45)
f_los_m = np.zeros(45)
f_los_med = np.zeros(45)
f_turb_m = np.zeros(45)
f_turb_med = np.zeros(45)
b_temp_m = np.zeros(45)
b_temp_med = np.zeros(45)
b_los_m = np.zeros(45)
b_turb_m = np.zeros(45)
b_los_med = np.zeros(45)
b_turb_med = np.zeros(45)

#chi2
pxnl = np.min([f_n])
chi2_f = np.zeros(pxnl)

for px in range(0,pxnl):
    chi2_f[px] = np.sum(((f_o.dat[0,0,px,:,ss] - f_i.dat[0,0,px,:,ss])/f_o.weights[:,ss])**2)/len(f_wav)
chi2_fm = np.mean(chi2_f)

# maps
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
fn = 28 #frame no.
#file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
#file_ca8 = file_search(datadir,'crispex*'+pref[1]+'*.fcube')
file_cak = file_search(datadir,'crispex*'+pref[2]+'*.fcube')
#file_ha = file_search(datadir,'crispex*'+pref[3]+'*.fcube')
#cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])
#cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])
cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])
#cube_ha = lp_read(datadir+file_ha[0],datadir+file_ha[1])
c_map = cube_cak[fn,0,-1,:,:]
cak_int_un = (lp_read_scan(datadir+'cak_int_un.fcube'))[fn,:,:]

#wav = s.air2vac(i.wav)
if(ss == 0):
    ymin = 0.05
    ymax = np.max(np.mean(f_obs, axis = 0))
elif(ss == 1 or ss == 2):
    ymin = -0.055
    ymax = 0.005
else:
    ymin = -0.03
    ymax = 0.03


################
# Variations along the fibril
################

# clipping the depth between log = [0,-7]
dep_init = 13 # -> logtau = -7
dep_fin = 132 # -> logtau = 0
f_depRF_clipped = f_depRF[dep_init: dep_fin+1]
depth_grid_sp = 17

# specify the array for the RF: fb/bg
# fibril 
temperature = np.transpose(f_m.temp[0,0,:,dep_init:dep_fin+1].squeeze()*1e-3)
vlos = np.transpose(f_m.vlos[0,0,:,dep_init:dep_fin+1].squeeze()*1e-5)
vmt = np.transpose(f_m.vturb[0,0,:,dep_init:dep_fin+1].squeeze()*1e-5)
# RFs
rf_t = f_rf.rf[0,0,:,0,:,:,0]
rf_vlos = f_rf.rf[0,0,:,1,:,:,0]
rf_vmt = f_rf.rf[0,0,:,2,:,:,0]
outname = outdir+'inv_pro_'+index+'_'+mode+'_cycle'+cyc+'length_withbg.pdf'

# bg
temperature_bg = np.transpose(b_m.temp[0,0,:,dep_init:dep_fin+1].squeeze()*1e-3)
vlos_bg = np.transpose(b_m.vlos[0,0,:,dep_init:dep_fin+1].squeeze()*1e-5)
vmt_bg = np.transpose(b_m.vturb[0,0,:,dep_init:dep_fin+1].squeeze()*1e-5)
# RFs
rf_t_bg = b_rf.rf[0,0,:,0,:,:,0]
rf_vlos_bg = b_rf.rf[0,0,:,1,:,:,0]
rf_vmt_bg = b_rf.rf[0,0,:,2,:,:,0]

# graphics
plt.close('all')

f = plt.figure(figsize = [5.70*1.7,6-0.16])
gs = gridspec.GridSpec(3, 2) # grid scale of the 1st col.
aspect = float(f_n_valid)/float(len(f_depRF_clipped))/2.5 #len(f_depRF_clipped)/f_n_valid # panels aspect ratio (y/x)
left, right, wspace, bottom, top, hspace = 0.05, 0.93, 0.02, 0.08, 0.95, 0.0
gs.update(left=left,
          right=right,
          wspace=wspace,
          bottom=bottom,
          top=top,
          hspace = hspace)

#======================
#                     FFIBRIL
#======================

ax1 = plt.subplot(gs[0,0],adjustable = 'box',aspect = 'equal')
ax2 = plt.subplot(gs[1,0],adjustable = 'box',aspect = 'equal')
ax3 = plt.subplot(gs[2,0],adjustable = 'box',aspect = 'equal')

# RF depth grids
ytick_pos = np.linspace(1,6,6,dtype='uint')*depth_grid_sp
ytick_lab = np.round(f_depRF_clipped[ytick_pos])
        
# temperature
#norm_t = np.max(np.abs(temperature))
temp_min, temp_max = 3.75, 9
panel1 = ax1.imshow(temperature, cmap = 'hot', vmin = temp_min, vmax = temp_max, aspect = aspect)

# max formation height
height_cak = np.zeros(f_n_valid,dtype=int)
height_ca8 = np.zeros(f_n_valid,dtype=int)
height_fe = np.zeros(f_n_valid,dtype=int)
for hh in range(f_n_valid):
    height_cak[hh] = np.amin(np.argwhere(rf_t[hh,:,0:143]>0)[:,0])
    height_ca8[hh] = np.amin(np.argwhere(rf_t[hh,:,146:234]>0)[:,0])
    height_fe[hh] = np.amin(np.argwhere(rf_t[hh,:,236:-1]>0)[:,0])
# plot max formation
ax12 = ax1.twinx()
ax12.set_xlim(0, temperature.shape[1]-1)
ax12.set_ylim(dep_fin, dep_init)
#ax12.plot(height_ca8, color = 'dimgray', linewidth = 1., label = r'FH$_{\rm8542, max}$')
#ax12.plot(height_cak, color = 'darkgray', linewidth = 1., label = r'FH$_{\rmcak, max}$')
#ax12.plot(dep_init+height_fe, color = 'green')
# plot RF max depth
RFmax_t = np.zeros(f_n_valid,dtype = int)
#print 'max Cak formation height = ', np.mean(f_dep[height_cak])
for pp in range(f_n_valid):
    # the most lower layers (that has response to cak)
    high_indx = 66
    RFmax_t[pp] = np.unravel_index(np.argmax(rf_t[pp,0:high_indx,:]), dims = rf_t[pp,0:high_indx,:].shape)[0]
ax12.plot(RFmax_t, color = 'black', linewidth = 1., label = r'RF$_{\rmcak, max}$')

ax12.yaxis.set_visible(False)

#f.legend(loc = 'upper center',
#                  fontsize = 8.5, handlelength = 0.9,framealpha = 0.5, ncol = 3,
 #                 frameon = False,
  #                #fancybox = False,
#)
    
#ax1.set_xlabel(r'Wavelength [$\mathrm{\AA}$]', fontdict=font)
ax1.set_ylabel(r'log($\tau_{500}$)')
ax1.set_xticklabels([])
ax1.set_yticks(ytick_pos)
ax1.set_yticklabels(ytick_lab.astype(int))
ax1.set_xlim(0, temperature.shape[1]-1)
for pp in range(gg):
    # pixel  indicator
    ax1.axvline(tp[pp], color = 'black', linestyle = '--', linewidth = 1.)
    ax1.axvline(tp[pp], color = 'white', linestyle = '--', linewidth = 1.)
    ax2.axvline(tp[pp], linestyle = '--', color = 'black', linewidth = 1.)
    ax3.axvline(tp[pp], linestyle = '--', color = 'white', linewidth = 1.)
    # pixel marker
    ax1.plot(tp[pp],6.5, marker = markers[pp], markersize = 6, color = 'red')
        

# vlos
norm_v = np.max(np.abs(vlos))
panel2 = ax2.imshow(vlos, cmap = 'bwr', vmin = -norm_v, vmax = norm_v,aspect = aspect)
ax2.set_ylabel(r'log($\tau_{500}$)')
ax2.set_xticklabels([])
ax2.set_yticks(ytick_pos)
ax2.set_yticklabels(ytick_lab.astype(int))

# vmt
norm_vmt = np.max(np.abs(vmt))
panel3 = ax3.imshow(vmt, cmap = 'bone', vmin = 0, vmax = norm_vmt,aspect = aspect)
ax3.set_xlabel(r'pixel along the fibril')
ax3.set_ylabel(r'log($\tau_{500}$)')
ax3.set_yticks(ytick_pos)
ax3.set_yticklabels(ytick_lab.astype(int))

plt.tight_layout()
plt.show()

#=========================
#                          Background
#=========================

ax1bg = plt.subplot(gs[0,1],adjustable = 'box',aspect = 'equal')
ax2bg = plt.subplot(gs[1,1],adjustable = 'box',aspect = 'equal')
ax3bg = plt.subplot(gs[2,1],adjustable = 'box',aspect = 'equal')

# RF depth grids
ytick_pos = np.linspace(1,6,6,dtype='uint')*depth_grid_sp
ytick_lab = np.round(f_depRF_clipped[ytick_pos])
        
# temperature_bg
#norm_t = np.max(np.abs(temperature_bg))
temp_min, temp_max = 3.75, 9
panel1 = ax1bg.imshow(temperature_bg, cmap = 'hot', vmin = temp_min, vmax = temp_max, aspect = aspect)

# max formation height
height_cak = np.zeros(f_n_valid,dtype=int)
height_ca8 = np.zeros(f_n_valid,dtype=int)
height_fe = np.zeros(f_n_valid,dtype=int)
for hh in range(f_n_valid):
    height_cak[hh] = np.amin(np.argwhere(rf_t_bg[hh,:,0:143]>0)[:,0])
    height_ca8[hh] = np.amin(np.argwhere(rf_t_bg[hh,:,146:234]>0)[:,0])
    height_fe[hh] = np.amin(np.argwhere(rf_t_bg[hh,:,236:-1]>0)[:,0])
# plot max formation
ax1bg2 = ax1bg.twinx()
ax1bg2.set_xlim(0, temperature_bg.shape[1]-1)
ax1bg2.set_ylim(dep_fin, dep_init)
#ax1bg2.plot(height_ca8, color = 'dimgray', linewidth = 1., label = r'FH$_{\rm8542, max}$')
#ax1bg2.plot(height_cak, color = 'darkgray', linewidth = 1., label = r'FH$_{\rmcak, max}$')
#ax1bg2.plot(dep_init+height_fe, color = 'green')
# plot RF max depth
RFmax_t = np.zeros(f_n_valid,dtype = int)
#print 'max Cak formation height = ', np.mean(f_dep[height_cak])
for pp in range(f_n_valid):
    # the most lower layers (that has response to cak)
    high_indx = 66
    RFmax_t[pp] = np.unravel_index(np.argmax(rf_t_bg[pp,0:high_indx,:]), dims = rf_t_bg[pp,0:high_indx,:].shape)[0]
ax1bg2.plot(RFmax_t, color = 'black', linewidth = 1., label = r'RF$_{\rmcak, max}$')

ax1bg2.yaxis.set_visible(False)

#ax1bg.set_xlabel(r'Wavelength [$\mathrm{\AA}$]', fontdict=font)
ax1bg.set_ylabel('')
ax1bg.set_xticklabels([])
ax1bg.set_yticks(ytick_pos)
ax1bg.set_yticklabels('')
ax1bg.set_xlim(0, temperature_bg.shape[1]-1)
for pp in range(gg):
    # pixel  indicator
    ax1bg.axvline(tp[pp], color = 'black', linestyle = '--', linewidth = 1.)
    ax1bg.axvline(tp[pp], color = 'white', linestyle = '--', linewidth = 1.)
    ax2bg.axvline(tp[pp], linestyle = '--', color = 'black', linewidth = 1.)
    ax3bg.axvline(tp[pp], linestyle = '--', color = 'white', linewidth = 1.)
    # pixel marker
    ax1bg.plot(tp[pp],6.5, marker = markers[pp], markersize = 6, color = 'gray')
        
# temp colorbar
axins = inset_axes(ax1bg, #axis
                   width="2%",  # width = 10% of parent_bbox width
                   height="90%",  # height : 50%
                   loc='center left',
                   bbox_to_anchor=(1.01, 0., 1, 1),
                   bbox_transform=ax1bg.transAxes,
                   borderpad=0,
)
ct_ticks = np.linspace(np.round(temp_min),np.round(temp_max),np.abs(np.round(temp_min)-np.round(temp_max))+1,dtype='uint')
ct = plt.colorbar(panel1, cax=axins, orientation="vertical",
                  ticks = ct_ticks,
)
ct.set_label(r'$T$ [kK]')
ct.ax.set_yticklabels(map(str,ct_ticks[:-1])+['>'+str(ct_ticks[-1])])

# vlos_bg
norm_v = np.max(np.abs(vlos_bg))
panel2 = ax2bg.imshow(vlos_bg, cmap = 'bwr', vmin = -norm_v, vmax = norm_v,aspect = aspect)
ax2bg.set_ylabel('')
ax2bg.set_xticklabels([])
ax2bg.set_yticks(ytick_pos)
ax2bg.set_yticklabels('')

# temp colorbar
axins = inset_axes(ax2bg, #axis
                   width="2%",  # width = 10% of parent_bbox width
                   height="90%",  # height : 50%
                   loc='center left',
                   bbox_to_anchor=(1.01, 0., 1, 1),
                   bbox_transform=ax2bg.transAxes,
                   borderpad=0,
)
cv = plt.colorbar(panel2, cax=axins, orientation="vertical")
cv.set_label(r'$v_{\rm LOS}$ [km s$^{-1}$]')

# vmt_bg
norm_vmt_bg = np.max(np.abs(vmt_bg))
panel3 = ax3bg.imshow(vmt_bg, cmap = 'bone', vmin = 0, vmax = norm_vmt_bg,aspect = aspect)
ax3bg.set_xlabel(r'pixel along the background')
ax3bg.set_ylabel('')
ax3bg.set_yticks(ytick_pos)
ax3bg.set_yticklabels('')

# vmt_bg colorbar
# temp colorbar
axins = inset_axes(ax3bg, #axis
                   width="2%",  # width = 10% of parent_bbox width
                   height="90%",  # height : 50%
                   loc='center left',
                   bbox_to_anchor=(1.01, 0., 1, 1),
                   bbox_transform=ax3bg.transAxes,
                   borderpad=0,
)
cvmt_bg = plt.colorbar(panel3, cax=axins, orientation="vertical")
cvmt_bg.set_label(r'$v_{\rm turb}$ [km s$^{-1}$]')
plt.tight_layout()
plt.show()

print '=============='
#f.savefig(outname, quality = 100)
print 'file saved to:'+outname
