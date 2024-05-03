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
f_rf = sp.profile('/home/seki2695/INV/stic/fb_back/f'+mode+'_synthetic_cycle2_'+index+'.nc')
# backgraound
b_i = sp.profile('/home/seki2695/INV/stic/bg_back/b'+mode+'_observed_'+index+'.nc')
b_o = sp.profile('/home/seki2695/INV/stic/bg_back/b'+mode+'_synthetic_cycle'+cyc+'_'+index+'.nc')
b_m = sp.model('/home/seki2695/INV/stic/bg_back/b'+mode+'_atmosout_cycle'+cyc+'_'+index+'.nc')
b_rf = sp.profile('/home/seki2695/INV/stic/bg_back/b'+mode+'_synthetic_cycle2_'+index+'.nc')
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


# Plotting RF #
##########

# fibril RF
rf_t = f_rf.rf[0,0,:,0,:,0:143,0]
rf_vlos = f_rf.rf[0,0,:,1,:,0:143,0]
rf_vmt = f_rf.rf[0,0,:,2,:,0:143,0]
outname = outdir+'inv_RF_'+index+'_'+mode+'_cycle'+cyc+'modified.pdf'
temp_color_map = 'afmhot'
marker_color = 'red'

# setting the graphics window
plt.close('all')
f = plt.figure(figsize = [3*gg,3*gg])

left = 0.055
right = 0.93
wspace = 0
bottom = 0.05
top = 0.95
hspace = 0

gs = gridspec.GridSpec(3, 3) # grid scale of the 2nd & 3rd col.
gs.update(left=left,
          right=right,
          wspace = wspace,
          bottom=bottom,
          top=top,
          hspace=hspace)

for pp in range(gg):
        
    # Vlos RF
    ax2 = host_subplot(gs[pp,0], adjustable='box', aspect='equal')
    norm_vlos = np.max(np.abs(rf_vlos[tp[pp],:,:]))
    vlos = ax2.imshow(im.histo_opt(rf_vlos[tp[pp],:,:]),cmap = 'RdGy',aspect = 'equal',vmin = -norm_vlos, vmax = norm_vlos)
    
    # pixel marker
    ax2.plot(135,140, marker = markers[pp], markersize = 7, color = marker_color, markeredgecolor = 'black', markeredgewidth = 0.5)
    
    ax2.set_xlabel(r'$\lambda$ [$\mathrm{\AA}$]', fontdict=font)
    ax2.set_ylabel(r'log($\tau\mathrm{_{500}}$)', fontdict=font)
    xtick_pos = [14,72,130]
    xtick_lab = np.round(f_wav[xtick_pos], decimals = 1)
    ytick_pos = [13, 30, 47, 64, 81, 98, 115, 132]#[8,18,28,38,48,58,68,78,88]
    ytick_lab = np.round(f_depRF[ytick_pos])
    ax2.set_yticks(ytick_pos)
    ax2.set_yticklabels(ytick_lab.astype(int), fontdict=font)
    
    if pp==2:
        ax2.set_xticks(xtick_pos)
        ax2.set_xticklabels(xtick_lab, fontdict=font)
    else:
        ax2.set_xticklabels('')
        ax2.set_xlabel('')
        ax2.set_xticks([])
        
        
    # extra axis for vlos prof
    ax23 = ax2.twiny()
    f_vlosRF =f_mRF.vlos[0,0,tp[pp],:]
    b_vlosRF =b_mRF.vlos[0,0,tp[pp],:]
    ax23y = np.linspace(0,len(f_depRF)-1,147)
    ax23.set_xlim(-7.4,7.4)
    ax23.plot(f_vlosRF.squeeze()*1e-5,ax23y, 'k-',color = 'red',linewidth = 1.5, label = r'v$_{\mathrm{f}}$')
    ax23.plot(b_vlosRF.squeeze()*1e-5,ax23y, 'k-',color = 'gray',linewidth = 1.5, label = r'v$_{\mathrm{b}}$')
    
    if pp==0:
        ax23.set_xlabel(r'$v \mathrm{_{LOS}}$ [km/s]', fontdict=font)
    else:
        ax23.set_xlabel('')
        ax23.set_xticklabels('')
        ax23.set_xticks([])
        
    #Required to remove some white border               
    ax2.autoscale(axis = 'both', enable = True)
    
    # Vmt RF
    ax3 = host_subplot(gs[pp,1], adjustable='box', aspect='equal')
    zero_vmt = np.abs(np.min(rf_vmt[tp[pp],:,:]))
    norm_vmt = np.max(np.abs(rf_vmt[tp[pp],:,:]))# + zero_vmt
    vmt = ax3.imshow(im.histo_opt(rf_vmt[tp[pp],:,:]),cmap = 'RdGy',aspect = 'equal', vmin = -norm_vmt,  vmax = norm_vmt)
        
    # pixel marker
    ax3.plot(135,140, marker = markers[pp], markersize = 7, color = marker_color, markeredgecolor = 'black', markeredgewidth = 0.5)
    
    ax3.set_xlabel(r'$\lambda$ [$\mathrm{\AA}$]', fontdict=font)
    ax3.set_ylabel('')
    ax3.set_yticks(ytick_pos)
    ax3.set_yticklabels('')
    xtick_pos = [14,72,130]
    xtick_lab = np.round(f_wav[xtick_pos], decimals = 1)
    ax3.set_xticks(xtick_pos)
    ax3.set_xticklabels(xtick_lab, fontdict=font)
    ytick_pos = [13, 30, 47, 64, 81, 98, 115, 132]#[8,18,28,38,48,58,68,78,88]
    ytick_lab = np.round(f_depRF[ytick_pos])
    #ax3.set_yticks(ytick_pos)
    #ax3.set_yticklabels(ytick_lab.astype(int), fontdict=font)
    
    if pp==2:
        ax3.set_xticks(xtick_pos)
        ax3.set_xticklabels(xtick_lab, fontdict=font)
    else:
        ax3.set_xticklabels('')
        ax3.set_xlabel('')
        ax3.set_xticks([])
        
        
    # extra axis for vturb prof
    ax33 = ax3.twiny()
    f_vmtRF =f_mRF.vturb[0,0,tp[pp],:]
    b_vmtRF =b_mRF.vturb[0,0,tp[pp],:]
    ax33y = np.linspace(0,len(f_depRF)-1,147)
    ax33.set_xlim(-0.05,8)
    ax33.plot(f_vmtRF.squeeze()*1e-5,ax33y, 'k-',color = 'red',linewidth = 1.5, label = r'v$_{\rm f}$')
    ax33.plot(b_vmtRF.squeeze()*1e-5,ax33y, 'k-',color = 'gray',linewidth = 1.5, label = r'v$_{\rm b}$')
        
    if pp==0:
        ax33.set_xlabel(r'$v\mathrm{_{turb}}$ [km/s]', fontdict=font)
    else:
        ax33.set_xlabel('')
        ax33.set_xticklabels('')
        ax33.set_xticks([])
        
        
    #Required to remove some white border               
    ax3.autoscale(axis = 'both', enable = True)

    # Temperature RF
    ax1 = host_subplot(gs[pp, 2], adjustable='box', aspect='equal')
    zero_t = np.abs(np.min(rf_t[tp[pp],:,:]))
    norm_t = np.max(np.abs(rf_t[tp[pp],:,:]))#+zero_t
    temperature = ax1.imshow(im.histo_opt(rf_t[tp[pp],:,:]),cmap = 'Greys', vmin = 0, vmax = norm_t)

    high_indx = 66
    rf_max_coords = np.unravel_index(np.argmax(rf_t[tp[pp],0:high_indx,:]), dims = rf_t[tp[pp],0:high_indx,:].shape)
    ax1.axhline(rf_max_coords[0], color = 'gray', alpha = 0.25, linewidth = 3.)
    ax1.axvline(rf_max_coords[1], color = 'gray', alpha = 0.25, linewidth = 3.)

    # pixel marker
    ax1.plot(135,140, marker = markers[pp], markersize = 7., color = marker_color, markeredgecolor = 'black', markeredgewidth = 0.5)
    
    # RF max
    #ax1.plot(np.argmax(rf_t[tp[pp],:,:], axis = 0), color = 'white', linewidth = 0.8, linestyle = '--', label = r'R$_{\rm\lambda,max}$')
    # the most lower layers (that has response to cak)
    high_indx = 66
    rf_max_coords = np.unravel_index(np.argmax(rf_t[tp[pp],0:high_indx,:]), dims = rf_t[tp[pp],0:high_indx,:].shape)
    print 'RFmax depth = ', f_depRF[rf_max_coords[0]]
    
    ax1.set_ylabel('')
    xtick_pos = [14,72,130]
    xtick_lab = np.round(f_wav[xtick_pos], decimals = 1)
    ax1.set_xticks(xtick_pos)
    ax1.tick_params(axis='both', which='major', labelsize = 8)
    
    if pp==2:
        ax1.set_xticklabels(xtick_lab, fontdict=font)
        ax1.set_xlabel(r'$\lambda$ [$\mathrm{\AA}$]', fontdict=font)
        #ax1.spines['top'].set_color('white')
    else:
        ax1.set_xticklabels('')
        ax1.set_xlabel('')
        ax1.set_xticks([])
        #if pp==1:
            #ax1.spines['top'].set_color('white')
        #ax1.spines['bottom'].set_color('white')        


    ytick_pos = [13, 30, 47, 64, 81, 98, 115, 132]#[8,18,28,38,48,58,68,78,88]
    ytick_lab = np.round(f_depRF[ytick_pos])
    ax1.set_yticks(ytick_pos)
    ax1.set_yticklabels('')
    
    # extra axis for intensity prof
    ax12 = ax1.twinx()
    ax12.spines['bottom'].set_color('white')
    ax12.set_ylim(0.01, 0.399)
    ax12.set_ylabel(r'$I$ [erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$]', fontdict=font)
    # BG intensity profile
    ax12.plot(b_obs[tp[pp],0:142],'.',color = 'darkgray', marker = '.',markersize = 3.5)#, label = r'I$_{\rm obs}$')
    # BG observed points
    ax12.plot(b_syn[tp[pp],0:142], color='darkgray',linewidth = 0.75, label = r'I$_{\rmb}$', linestyle = '--')
    # FB intensity profile
    ax12.plot(f_obs[tp[pp],0:142],'.',color = 'red', marker = '.',markersize = 3.5)#, label = r'I$_{\rm obs}$')
    # FB observed points
    ax12.plot(f_syn[tp[pp],0:142], color='r',linewidth = 0.75, label = r'I$_{\rmf}$', linestyle = '--')
    
    # extra axis for temperature prof
    ax13 = ax1.twiny()
    f_tempRF =f_mRF.temp[0,0,tp[pp],:]
    b_tempRF =b_mRF.temp[0,0,tp[pp],:]
    ax13y = np.linspace(0,len(f_depRF)-1,147)
    ax13.set_xlim(3.75,9)
    
    if pp==0:
        ax13.set_xlabel(r'$T$ [kK]', fontdict=font)
    else:
        ax13.set_xlabel('')
        ax13.set_xticklabels('')
        ax13.set_xticks([])
            

        
    # bg temperature profile
    ax13.plot(b_tempRF.squeeze()*1e-3,ax13y, 'k-',color = 'gray',linewidth = 1.5, label = r'T$_{\rmb}$')
    # fibril temperture profile
    ax13.plot(f_tempRF.squeeze()*1e-3,ax13y, 'k-',color = 'red',linewidth = 1.5, label = r'T$_{\rmf}$')
    
    #Required to remove some white border               
    ax1.autoscale(axis = 'both', enable = True)

    ax1.tick_params(axis='both', which='major', labelsize = 8)
    ax2.tick_params(axis='both', which='major', labelsize = 8)
    ax3.tick_params(axis='both', which='major', labelsize = 8)
    ax12.tick_params(axis='both', which='major', labelsize = 8)
    ax13.tick_params(axis='both', which='major', labelsize = 8)
    ax23.tick_params(axis='both', which='major', labelsize = 8)
    ax33.tick_params(axis='both', which='major', labelsize = 8)

    
    # legend settings
    if pp==0:
        leg = ax1.legend( loc='upper center',
                          ncol=5, borderaxespad=0., fontsize = 8.5,
                          frameon = False, #framealpha = 0.5,#fancybox = False,
                          handlelength = 0.9,
                          handletextpad = 0.1,
        )
        for text in leg.get_texts():
            plt.setp(text, color = 'black')
        
            ax2.legend(ncol = 2, loc = 'upper left', fontsize = 8.5, borderaxespad=0,
                       frameon = False, fancybox = False,
                       handlelength = 0.9,
                       handletextpad = 0.3,
                       columnspacing = 0.7,
            )
            ax3.legend(ncol = 2, loc = 'upper left', fontsize = 8.5, borderaxespad=0,
                       frameon = False, fancybox = False,
                       handlelength = 0.7,
                       handletextpad = 0.2,
                       columnspacing = 0.6,
            )

    if pp == 2:
        # temperature colorbar
        axins1 = inset_axes(ax1,
                            width="19%",  # width = 90% of parent_bbox width
                            height="5%",  # height : 50%
                            loc='lower left')
        # change the color of the colomap ticks to white
        axins1.axes.tick_params(which = 'major', length = 3, color = 'black', #tick color
                                labelcolor = 'black', labelsize = 8) #tick label color and font size
        ct = plt.colorbar(temperature, cax=axins1, orientation="horizontal", ticks=[0,norm_t])
        ct.ax.xaxis.tick_top()
        ct.ax.set_xticklabels(['0','1'])
        ct.outline.set_edgecolor('black')
        
        # vlos colorbar
        axins2 = inset_axes(ax2,
                            width="38%",  # width = 90% of parent_bbox width
                            height="5%",  # height : 50%
                            loc='lower left')
        # change the color of the colomap ticks to white
        axins2.axes.tick_params(which = 'major', length = 3, labelsize = 8) #tick label color and font size
        cv = plt.colorbar(vlos, cax=axins2, orientation="horizontal", ticks=[-norm_vlos,0,norm_vlos])
        cv.ax.xaxis.tick_top()
        cv.ax.set_xticklabels(['-1','0','1'])
            
        # vmt colorbar
        axins3 = inset_axes(ax3,
                            width="38%",  # width = 90% of parent_bbox width
                            height="5%",  # height : 50%
                            loc='lower left')
        # change the color of the colomap ticks to white
        axins3.axes.tick_params(which = 'major', length = 3, labelsize = 8) #tick label color and font size
        cvmt = plt.colorbar(vmt, cax=axins3, orientation="horizontal", ticks=[-norm_vmt,0,norm_vmt])
        cvmt.ax.xaxis.tick_top()
        cvmt.ax.set_xticklabels(['-1','0','1'])
        
# f.tight_layout()
plt.show()
#f.savefig(outname, quality = 100)
print 'file saved to: '+outname

