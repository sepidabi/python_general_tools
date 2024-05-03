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


# Plotting the fibril profile #
#################
if plot_mode == 1:
    plt.close("all")

    f = plt.figure(figsize=(8,8))
    ax1 = plt.subplot2grid((10,6), (0,0), colspan=6, rowspan=3)
    ax2 = plt.subplot2grid((10,6), (3,0), colspan=2, rowspan=3)
    ax3 = plt.subplot2grid((10,6), (3,2), colspan=2, rowspan=3)
    ax4 = plt.subplot2grid((10,6), (3,4), colspan=2, rowspan=3)
    ax5 = plt.subplot2grid((10,6),(6,0), colspan=3,rowspan=6)
    ax6 = plt.subplot2grid((10,6),(6,3), colspan=3,rowspan=6)

    # amplifying the cak spectral profile
    amp = 3.5 #cak amplifier
    f_obs_amp = f_i.dat[0,0,:,:,ss]
    f_obs_amp[:,0:42]  = f_i.dat[0,0,:,0:42,ss]*amp
    f_syn_amp = f_o.dat[0,0,:,:,ss]
    f_syn_amp[:,0:42]  = f_o.dat[0,0,:,0:42,ss]*amp
    b_obs_amp = b_i.dat[0,0,:,:,ss]
    b_obs_amp[:,0:42]  = b_i.dat[0,0,:,0:42,ss]*amp
    b_syn_amp = b_o.dat[0,0,:,:,ss]
    b_syn_amp[:,0:42]  = b_o.dat[0,0,:,0:42,ss]*amp


    ax1.plot(np.mean(b_obs,axis = 0),'s',color = 'blue',markersize = 3)
    ax1.plot(np.mean(f_obs, axis = 0),'.',color = 'red',markersize = 4.5)
    ax1.plot(np.mean(b_syn,axis = 0), color='b',linewidth = 1.25)
    ax1.plot(np.mean(f_syn, axis = 0), color='r',linewidth = 1.25)
    
    for ff in range(0,f_n):
        f_temp_m += f_m.temp[0,0,ff,:]
        ax2.plot(f_dep, f_m.temp[0,0,ff,:].squeeze()*1e-3, 'k-',color = 'red',alpha = alp)
        ax3.plot(f_dep, f_m.vturb[0,0,ff,:].squeeze()*1.e-5, 'k-',color = 'red',alpha = alp)
        ax4.plot(f_dep, f_m.vlos[0,0,ff,:].squeeze()*1.e-5, 'k-',color = 'red',alpha = alp)
    for dd in range(45):
        f_temp_m[dd] = np.mean(f_m.temp[0,0,:,dd])
        f_temp_med[dd] = np.median(f_m.temp[0,0,:,dd])
        f_los_m[dd] = np.mean(f_m.vlos[0,0,:,dd])
        f_los_med[dd] = np.median(f_m.vlos[0,0,:,dd])
        f_turb_m[dd] = np.mean(f_m.vturb[0,0,:,dd])
        f_turb_med[dd] = np.median(f_m.vturb[0,0,:,dd])

    for bb in range(0,b_n):
        b_temp_m += b_m.temp[0,0,bb,:]
        ax2.plot(b_dep, b_m.temp[0,0,bb,:].squeeze()*1e-3, 'k-',color = 'b',alpha = alp)
        ax3.plot(b_dep, b_m.vturb[0,0,bb,:].squeeze()*1.e-5, 'k-',color = 'b',alpha = alp)
        ax4.plot(b_dep, b_m.vlos[0,0,bb,:].squeeze()*1.e-5, 'k-',color = 'b',alpha = alp)
    for dd in range(45):
        b_temp_m[dd] = np.mean(b_m.temp[0,0,:,dd])
        b_temp_med[dd] = np.median(b_m.temp[0,0,:,dd])
        b_los_m[dd] = np.mean(b_m.vlos[0,0,:,dd])
        b_los_med[dd] = np.median(b_m.vlos[0,0,:,dd])
        b_turb_m[dd] = np.mean(b_m.vturb[0,0,:,dd])
        b_turb_med[dd] = np.median(b_m.vturb[0,0,:,dd])
 
    
    lw = 1.1 #plot linewidth
    ax2.plot(b_dep, (b_temp_m).squeeze()*1e-3, 'k-',color = 'b',linewidth = lw)
    ax2.plot(b_dep, (b_temp_med).squeeze()*1e-3, 'k-',color = 'b',linestyle = '--',linewidth = lw)
    ax2.plot(f_dep, (f_temp_m).squeeze()*1e-3, 'k-',color = 'r',linewidth = lw)
    ax2.plot(f_dep, (f_temp_med).squeeze()*1e-3, 'k-',color = 'r',linestyle = '--',linewidth = lw)
    for n in range(len(nodes_temp)):
        ax2.axvline(nodes_temp[n],linestyle = '--',color = 'k',linewidth = 0.2)
        ax3.plot(b_dep, (b_turb_m).squeeze()*1e-5, 'k-',color = 'b',linewidth = lw)
        ax3.plot(b_dep, (b_turb_med).squeeze()*1e-5, 'k-',color = 'b',linestyle = '--',linewidth = lw)
        ax3.plot(f_dep, (f_turb_m).squeeze()*1e-5, 'k-',color = 'r',linewidth = lw)
        ax3.plot(f_dep, (f_turb_med).squeeze()*1e-5, 'k-',color = 'r',linestyle = '--',linewidth = lw)
    for n in range(len(nodes_vturb)):
        ax3.axvline(nodes_vturb[n],linestyle = '--',color = 'k',linewidth = 0.2)
        ax4.plot(b_dep, (b_los_m).squeeze()*1e-5, 'k-',color = 'b',linewidth = lw)
        ax4.plot(b_dep, (b_los_med).squeeze()*1e-5, 'k-',color = 'b',linestyle = '--',linewidth = lw)
        ax4.plot(f_dep, (f_los_m).squeeze()*1e-5, 'k-',color = 'r',linewidth = lw)
        ax4.plot(f_dep, (f_los_med).squeeze()*1e-5, 'k-',color = 'r',linestyle = '--',linewidth = lw)
    for n in range(len(nodes_vlos)):
        ax4.axvline(nodes_vlos[n],linestyle = '--',color = 'k',linewidth = 0.2)

    ax1.set_ylabel(r'I$_{mean}$ [erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$]')
    ax1.set_xlabel(r'Wavelength [$\mathrm{\AA}$]')
    ax1.set_xticks(([21,86,153,211]))
    ax1.set_xticklabels((['3934','8542','6301','6302']))
    ax1.set_ylim(ymin, ymax+0.05)
    ax1.set_xlim(0,b_wav.size-1)
    ax1.legend(['Background','Fibril'])
    
    ax2.set_ylim(3.75,11)
    ax2.set_xlim(-6.5,1)
    ax2.set_xlabel(r'log($\mathrm{\tau_{500}}$)')
    ax2.set_ylabel(r'T [kK]')

    ax3.set_ylim(-0.1,10)
    ax3.set_xlim(-7.5,1)
    ax3.set_xlabel(r'log($\mathrm{\tau_{500}}$)')
    ax3.set_ylabel(r'v$_{\rm turb}$ [km/s]')

    ax4.set_ylim(-7.5,7.5)
    ax4.set_xlim(-7.5,1)
    ax4.set_xlabel(r'log($\mathrm{\tau_{500}}$)')
    ax4.set_ylabel(r'v$\mathrm{_{LOS}}$ [km/s]')

    #inv = restore('inv_in'+index+'.sav')
    #xmin = 0
    #xmax = c_map.shape[1]
    #ymin = 0
    #ymax = c_map.shape[0]
    xmin = int(np.min([f_x_pts[0],f_x_pts[-1]])-20)
    xmax = int(np.max([f_x_pts[0],f_x_pts[-1]])+20)
    ymin = int(np.min([f_y_pts[0],f_y_pts[-1]])-20)
    ymax = int(np.max([f_y_pts[0],f_y_pts[-1]])+20)
    xlabel = '[pixels]'
    ylabel = '[pixels]'
    cmap = 'gray'
    title = r'Continuum 4000 $\mathrm{\AA}$'
    titlek = r'Ca II K $\mathrm{\lambda}$-integrated'

    # cont map
    c_map_cropped = c_map[ymin:ymax,xmin:xmax]
    ax5.imshow(c_map_cropped, cmap = cmap)
    ax5.set_xlim(0, xmax-xmin)
    ax5.set_ylim(0, ymax-ymin)
    ax5.set_xlabel(xlabel)
    ax5.set_ylabel(ylabel)
    #ax5.set_title(title)
    ax5.text(0.5,0.5,r'Continuum 4000 $\mathrm{\AA}$',color = 'white',horizontalalignment='left', verticalalignment='bottom') #, bbox={'facecolor': 'black', 'pad': 10,'alpha':0.5}

    # vlos profile along the fibril
    test = np.zeros(c_map_cropped.shape)

    # overplotting the background
    ax5.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue')
    # overplotting the fibril
    ax5.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'red')
    # marking the extreme pixels along the fibril
    for ii in range(gg):
        ax5.plot(b_x_pts[tp[ii]]-xmin,b_y_pts[tp[ii]]-ymin, marker = markers[ii], color = 'blue',markersize = 5)
        ax5.plot(f_x_pts[tp[ii]]-xmin,f_y_pts[tp[ii]]-ymin, marker = markers[ii], color = 'red',markersize = 5)
        ax6.plot(b_x_pts[tp[ii]]-xmin,b_y_pts[tp[ii]]-ymin, marker = markers[ii], color = 'blue',markersize = 5,markerfacecolor = 'none')
        ax6.plot(f_x_pts[tp[ii]]-xmin,f_y_pts[tp[ii]]-ymin, marker = markers[ii], color = 'red',markersize = 5,markerfacecolor = 'none')
    # annotation board
    ax5.text(0.5,120,reg, color = 'black', fontsize = 9)

    # cak map
    cak_map_cropped = cak_int_un[ymin:ymax,xmin:xmax]
    ax6.imshow(cak_map_cropped, cmap = cmap)
    ax6.set_xlim(0, xmax-xmin)
    ax6.set_ylim(0, ymax-ymin)
    ax6.set_xlabel(xlabel)
    ax6.set_ylabel('')
    #ax6.set_title(titlek)
    ax6.text(0.5,0.5,r'Ca II K $\mathrm{\lambda}$-integrated',color = 'white',horizontalalignment='left', verticalalignment='bottom') #, bbox={'facecolor': 'black', 'pad': 10,'alpha':0.5}
    #ax6.text(7,95,'227',color = 'white',horizontalalignment='left', verticalalignment='bottom') #, bbox={'facecolor': 'black', 'pad': 10,'alpha':0.5}

    # overplotting the fibril
    ax6.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue',alpha = 0.3,linewidth = 4)
    ax6.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'red',alpha = 0.3,linewidth = 4)
    f.set_tight_layout(True)
    plt.ion()
    plt.show()

    f.savefig(outdir+'inv_pro_'+index+'_'+mode+'_cycle'+cyc+'.pdf', quality = 100)
    print 'file saved to:'+outdir+'inv_pro_'+index+'_'+mode+'_cycle'+cyc+'.pdf'

    #plt.close('all')
    
# Plotting RF #
##########
if plot_mode == 2:

    # fibril RF
    rf_t = f_rf.rf[0,0,:,0,:,0:143,0]
    rf_vlos = f_rf.rf[0,0,:,1,:,0:143,0]
    rf_vmt = f_rf.rf[0,0,:,2,:,0:143,0]
    outname = outdir+'inv_RF_'+index+'_'+mode+'_cycle'+cyc+'modified.pdf'
    temp_color_map = 'afmhot'
    marker_color = 'red'
    
    # setting the graphics window
    plt.close('all')
    f = plt.figure(figsize = [3*gg+0.3,3*gg-0.14])

    gs1 = gridspec.GridSpec(3, 1) # grid scale of the 1st col.
    gs1.update(left=0.045,
               right=0.318,
               wspace=0.,
               bottom=0.05,
               top=0.94,
               hspace = 0.)
    gs = gridspec.GridSpec(3, 2) # grid scale of the 2nd & 3rd col.
    gs.update(left=0.44,
              right=0.99,
              wspace = 0.,
              bottom=0.05,
              top=0.94,
              hspace=0.)
    
    for pp in range(gg):

        # Temperature RF
        ax1 = host_subplot(gs1[pp, 0], adjustable='box', aspect='equal')
        zero_t = np.abs(np.min(rf_t[tp[pp],:,:]))
        norm_t = np.max(np.abs(rf_t[tp[pp],:,:]))#+zero_t
        temperature = ax1.imshow(im.histo_opt(rf_t[tp[pp],:,:]),cmap = temp_color_map, vmin = 0, vmax = norm_t)
        # pixel marker
        ax1.plot(135,140, marker = markers[pp], markersize = 6, color = marker_color, markeredgecolor = 'white', markeredgewidth = 0.5)
        
        # RF max
        ax1.plot(np.argmax(rf_t[tp[pp],:,:], axis = 0), color = 'white', linewidth = 0.8, linestyle = '--', label = r'R$_{\rm\lambda,max}$')
        # the most lower layers (that has response to cak)
        high_indx = 66
        rf_max_coords = np.unravel_index(np.argmax(rf_t[tp[pp],0:high_indx,:]), dims = rf_t[tp[pp],0:high_indx,:].shape)
        print 'RFmax depth = ', f_depRF[rf_max_coords[0]]
        ax1.axhline(rf_max_coords[0], color = 'white', alpha = 0.25, linewidth = 3.)
        ax1.axvline(rf_max_coords[1], color = 'white', alpha = 0.25, linewidth = 3.)
        
        ax1.set_ylabel(r'log($\mathrm{\tau_{500}}$)', fontdict=font)
        xtick_pos = [14,72,130]
        xtick_lab = np.round(f_wav[xtick_pos], decimals = 1)
        ax1.set_xticks(xtick_pos)
        ax1.tick_params(axis='both', which='major', labelsize = 8)
        
        if pp==2:
            ax1.set_xticklabels(xtick_lab, fontdict=font)
            ax1.set_xlabel(r'Wavelength [$\mathrm{\AA}$]', fontdict=font)
            ax1.spines['top'].set_color('white')
        else:
            ax1.set_xticklabels('')
            ax1.set_xlabel('')
            ax1.set_xticks([])
            if pp==1:
                ax1.spines['top'].set_color('white')
            ax1.spines['bottom'].set_color('white')        


        ytick_pos = [13, 30, 47, 64, 81, 98, 115, 132]#[8,18,28,38,48,58,68,78,88]
        ytick_lab = np.round(f_depRF[ytick_pos])
        ax1.set_yticks(ytick_pos)
        ax1.set_yticklabels(ytick_lab.astype(int), fontdict=font)
            
        # extra axis for intensity prof
        ax12 = ax1.twinx()
        ax12.spines['bottom'].set_color('white')
        ax12.set_ylim(0.01, 0.399)
        ax12.set_ylabel(r'I [erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$]', fontdict=font)
        # BG intensity profile
        ax12.plot(b_obs[tp[pp],0:142],'.',color = 'darkgray', marker = '.',markersize = 2.)#, label = r'I$_{\rm obs}$')
        # BG observed points
        ax12.plot(b_syn[tp[pp],0:142], color='darkgray',linewidth = 1., label = r'I$_{\rmb}$')
        # FB intensity profile
        ax12.plot(f_obs[tp[pp],0:142],'.',color = 'red', marker = '.',markersize = 2.)#, label = r'I$_{\rm obs}$')
        # FB observed points
        ax12.plot(f_syn[tp[pp],0:142], color='r',linewidth = 1., label = r'I$_{\rmf}$')

        # extra axis for temperature prof
        ax13 = ax1.twiny()
        f_tempRF =f_mRF.temp[0,0,tp[pp],:]
        b_tempRF =b_mRF.temp[0,0,tp[pp],:]
        ax13y = np.linspace(0,len(f_depRF)-1,147)
        ax13.set_xlim(3.75,9)

        if pp==0:
            ax13.set_xlabel(r'T [kK]', fontdict=font)
        else:
            ax13.set_xlabel('')
            ax13.set_xticklabels('')
            ax13.set_xticks([])
            

        
        # bg temperature profile
        ax13.plot(b_tempRF.squeeze()*1e-3,ax13y, 'k-',color = 'dimgray',linewidth = 1., label = r'T$_{\rmb}$')
        # fibril temperture profile
        ax13.plot(f_tempRF.squeeze()*1e-3,ax13y, 'k-',color = 'orange',linewidth = 1., label = r'T$_{\rmf}$')
        
        #Required to remove some white border               
        ax1.autoscale(axis = 'both', enable = True)
        
        # Vlos RF
        ax2 = host_subplot(gs[pp,0], adjustable='box', aspect='equal')
        norm_vlos = np.max(np.abs(rf_vlos[tp[pp],:,:]))
        vlos = ax2.imshow(im.histo_opt(rf_vlos[tp[pp],:,:]),cmap = 'RdGy',aspect = 'equal',vmin = -norm_vlos, vmax = norm_vlos)
        
        # pixel marker
        ax2.plot(135,140, marker = markers[pp], markersize = 6, color = marker_color, markeredgecolor = 'black', markeredgewidth = 0.5)
            
        ax2.set_xlabel(r'Wavelength [$\mathrm{\AA}$]', fontdict=font)
        ax2.set_ylabel(r'log($\mathrm{\tau_{500}}$)', fontdict=font)
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
        ax23.plot(f_vlosRF.squeeze()*1e-5,ax23y, 'k-',color = 'red',linewidth = 1., label = r'v$_{\mathrm{f}}$')
        ax23.plot(b_vlosRF.squeeze()*1e-5,ax23y, 'k-',color = 'gray',linewidth = 1., label = r'v$_{\mathrm{b}}$')

        if pp==0:
            ax23.set_xlabel(r'v$\mathrm{_{LOS}}$ [km/s]', fontdict=font)
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
        ax3.plot(135,140, marker = markers[pp], markersize = 6, color = marker_color, markeredgecolor = 'black', markeredgewidth = 0.5)
        
        ax3.set_xlabel(r'Wavelength [$\mathrm{\AA}$]', fontdict=font)
        ax3.set_ylabel('')
        ax3.set_yticks([])
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
        ax33.plot(f_vmtRF.squeeze()*1e-5,ax33y, 'k-',color = 'red',linewidth = 1., label = r'v$_{\rm f}$')
        ax33.plot(b_vmtRF.squeeze()*1e-5,ax33y, 'k-',color = 'gray',linewidth = 1., label = r'v$_{\rm b}$')

        ax1.tick_params(axis='both', which='major', labelsize = 8)
        ax2.tick_params(axis='both', which='major', labelsize = 8)
        ax3.tick_params(axis='both', which='major', labelsize = 8)
        ax12.tick_params(axis='both', which='major', labelsize = 8)
        ax13.tick_params(axis='both', which='major', labelsize = 8)
        ax23.tick_params(axis='both', which='major', labelsize = 8)
        ax33.tick_params(axis='both', which='major', labelsize = 8)

        if pp==0:
            ax33.set_xlabel(r'v$\mathrm{_{turb}}$ [km/s]', fontdict=font)
        else:
            ax33.set_xlabel('')
            ax33.set_xticklabels('')
            ax33.set_xticks([])

        
        #Required to remove some white border               
        ax3.autoscale(axis = 'both', enable = True)
        
        # legend settings
        if pp==0:
            leg = ax1.legend( loc='upper left',
                              ncol=5, borderaxespad=0., fontsize = 8.5,
                              frameon = False, #framealpha = 0.5,#fancybox = False,
                              handlelength = 0.7,
                              handletextpad = 0.1,
            )
            for text in leg.get_texts():
                plt.setp(text, color = 'w')
                    
            ax2.legend(ncol = 2, loc = 'upper left', fontsize = 8.5, borderaxespad=0,
                       frameon = False, fancybox = False,
                       handlelength = 0.7,
                       handletextpad = 0.2,
                       columnspacing = 0.6,
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
                                width="38%",  # width = 90% of parent_bbox width
                                height="5%",  # height : 50%
                                loc='lower left')
            # change the color of the colomap ticks to white
            axins1.axes.tick_params(which = 'major', length = 3, color = 'white', #tick color
                                    labelcolor = 'white', labelsize = 8) #tick label color and font size
            ct = plt.colorbar(temperature, cax=axins1, orientation="horizontal", ticks=[0,norm_t/2,norm_t])
            ct.ax.xaxis.tick_top()
            ct.ax.set_xticklabels(['0','0.5','1'])
            ct.outline.set_edgecolor('white')
        
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

################
# Variations along the fibril
################
if plot_mode ==3:

    # clipping the depth between log = [0,-7]
    dep_init = 13 # -> logtau = -7
    dep_fin = 132 # -> logtau = 0
    f_depRF_clipped = f_depRF[dep_init: dep_fin+1]
    depth_grid_sp = 17

    # specify the array for the RF: fb/bg
    for what in range(2):
        if what == 0:
            # fibril 
            temperature = np.transpose(f_m.temp[0,0,:,dep_init:dep_fin+1].squeeze()*1e-3)
            vlos = np.transpose(f_m.vlos[0,0,:,dep_init:dep_fin+1].squeeze()*1e-5)
            vmt = np.transpose(f_m.vturb[0,0,:,dep_init:dep_fin+1].squeeze()*1e-5)
            # RFs
            rf_t = f_rf.rf[0,0,:,0,:,:,0]
            rf_vlos = f_rf.rf[0,0,:,1,:,:,0]
            rf_vmt = f_rf.rf[0,0,:,2,:,:,0]
            outname = outdir+'inv_pro_'+index+'_'+mode+'_cycle'+cyc+'length.pdf'
        else:
            # bg
            temperature = np.transpose(b_m.temp[0,0,:,dep_init:dep_fin+1].squeeze()*1e-3)
            vlos = np.transpose(b_m.vlos[0,0,:,dep_init:dep_fin+1].squeeze()*1e-5)
            vmt = np.transpose(b_m.vturb[0,0,:,dep_init:dep_fin+1].squeeze()*1e-5)
            # RFs
            rf_t = b_rf.rf[0,0,:,0,:,:,0]
            rf_vlos = b_rf.rf[0,0,:,1,:,:,0]
            rf_vmt = b_rf.rf[0,0,:,2,:,:,0]
            outname = outdir+'inv_pro_'+index+'_'+mode+'_cycle'+cyc+'lengthbg.pdf'

        # graphics
        plt.close('all')

        f = plt.figure(figsize = [5.70,6])
        gs = gridspec.GridSpec(3, 1) # grid scale of the 1st col.
        aspect = float(f_n_valid)/float(len(f_depRF_clipped))/2.5 #len(f_depRF_clipped)/f_n_valid # panels aspect ratio (y/x)
        gs.update(left=0.085,
                  right=0.88,
                  #wspace=0.05,
                  bottom=0.08,
                  top=0.99,
                  hspace = 0.0)
        
        ax1 = plt.subplot(gs[0,0],adjustable = 'box',aspect = 'equal')
        ax2 = plt.subplot(gs[1,0],adjustable = 'box',aspect = 'equal')
        ax3 = plt.subplot(gs[2,0],adjustable = 'box',aspect = 'equal')
                
        # normal depth grids
        #ytick_pos = [4, 9, 14, 19, 24, 29, 34, 39, 44]#[8,18,28,38,48,58,68,78,88]
        #ytick_lab = np.round(f_dep[ytick_pos])
        
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
        ax12.plot(height_ca8, color = 'dimgray', linewidth = 1., label = r'FH$_{\rm8542, max}$')
        ax12.plot(height_cak, color = 'darkgray', linewidth = 1., label = r'FH$_{\rmcak, max}$')
        #ax12.plot(dep_init+height_fe, color = 'green')
        # plot RF max depth
        RFmax_t = np.zeros(f_n_valid,dtype = int)
        #print 'max Cak formation height = ', np.mean(f_dep[height_cak])
        for pp in range(f_n_valid):
            print 'pixel no:', pp
            print 'max Cak formation height at pixel = ', f_depRF[height_cak[pp]]
            # the most lower layers (that has response to cak)
            high_indx = 66
            RFmax_t[pp] = np.unravel_index(np.argmax(rf_t[pp,0:high_indx,:]), dims = rf_t[pp,0:high_indx,:].shape)[0]
            print 'max RF = ', f_depRF[RFmax_t[pp]]
            print '==============='
        ax12.plot(RFmax_t, color = 'white', linewidth = 1., label = r'RF$_{\rmcak, max}$')

        ax12.yaxis.set_visible(False)

        leg = ax12.legend(loc = 'lower center',
                          fontsize = 8.5, handlelength = 0.9,framealpha = 0.5, ncol = 3
                    #frameon = False,
                    #fancybox = False,
        )
        for text in leg.get_texts():
                    plt.setp(text, color = 'black')

        #ax1.set_xlabel(r'Wavelength [$\mathrm{\AA}$]', fontdict=font)
        ax1.set_ylabel(r'log($\mathrm{\tau_{500}}$)')
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
            ax1.plot(tp[pp],6.5, marker = markers[pp], markersize = 6, color = 'black')

        # temp colorbar
        axins = inset_axes(ax1, #axis
                           width="2%",  # width = 10% of parent_bbox width
                           height="90%",  # height : 50%
                           loc='center left',
                           bbox_to_anchor=(1.01, 0., 1, 1),
                           bbox_transform=ax1.transAxes,
                           borderpad=0,
        )
        ct_ticks = np.linspace(np.round(temp_min),np.round(temp_max),np.abs(np.round(temp_min)-np.round(temp_max))+1,dtype='uint')
        ct = plt.colorbar(panel1, cax=axins, orientation="vertical",
                          ticks = ct_ticks,
        )
        ct.set_label(r'T (kK)')
        ct.ax.set_yticklabels(map(str,ct_ticks[:-1])+['>'+str(ct_ticks[-1])])
        
        # vlos
        norm_v = np.max(np.abs(vlos))
        panel2 = ax2.imshow(vlos, cmap = 'bwr', vmin = -norm_v, vmax = norm_v,aspect = aspect)
        ax2.set_ylabel(r'log($\mathrm{\tau_{500}}$)')
        ax2.set_xticklabels([])
        ax2.set_yticks(ytick_pos)
        ax2.set_yticklabels(ytick_lab.astype(int))
        
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
        cv.set_label(r'v$_{\rm LOS}$ [km/s]')
        
        # vmt
        norm_vmt = np.max(np.abs(vmt))
        panel3 = ax3.imshow(vmt, cmap = 'bone', vmin = 0, vmax = norm_vmt,aspect = aspect)
        ax3.set_xlabel(r'pixel along the path')
        ax3.set_ylabel(r'log($\mathrm{\tau_{500}}$)')
        ax3.set_yticks(ytick_pos)
        ax3.set_yticklabels(ytick_lab.astype(int))
        
        # vmt colorbar
        # temp colorbar
        axins = inset_axes(ax3, #axis
                           width="2%",  # width = 10% of parent_bbox width
                           height="90%",  # height : 50%
                           loc='center left',
                           bbox_to_anchor=(1.01, 0., 1, 1),
                           bbox_transform=ax3.transAxes,
                           borderpad=0,
        )
        cvmt = plt.colorbar(panel3, cax=axins, orientation="vertical")
        cvmt.set_label(r'v$_{\rm turb}$ [km/s]')
        plt.tight_layout()
        plt.show()

        print '=============='
        #f.savefig(outname, quality = 100)
        print 'file saved to:'+outname

if plot_mode==4:

    plt.close('all')
    # max height formation of the FIBRIL and BG
    fmh, bmh = 10, 11
    # height of max response function
    fhm, bhm = 16, 16

    # plot temperature in above heights
    plt.plot(f_m.temp[0,0,:,fmh].squeeze()*1e-3, color = 'darkred', label = 'F maximum formation height')
    plt.plot(b_m.temp[0,0,:,bmh].squeeze()*1e-3, color = 'darkblue', label = 'B maximum formation height')
    plt.plot(f_m.temp[0,0,:,fhm].squeeze()*1e-3, color = 'red', label = 'F RFmax')
    plt.plot(b_m.temp[0,0,:,bhm].squeeze()*1e-3, color = 'blue', label = 'B RFmax')
    plt.legend()
    plt.xlabel('pixel along the path')
    plt.ylabel('T [kK]')
    plt.xlim([0,f_l])
    for pp in range(gg):
        # pixel  indicator
        plt.axvline(tp[pp], color = 'black', linestyle = '--', linewidth = 1.)
        #ax2.axvline(tp[pp], linestyle = '--', color = 'black', linewidth = 1.)
        #ax3.axvline(tp[pp], linestyle = '--', color = 'black', linewidth = 1.)
        # pixel marker
        plt.plot(tp[pp],4, marker = markers[pp], markersize = 6, color = 'black')

    plt.tight_layout()
    plt.show()
    #plt.savefig(outdir+'inv_pro_'+index+'_tem_dif.pdf', quality = 100)
    print 'file saved to: ' + outdir+'inv_pro_'+index+'_tem_dif.pdf'

#plt.close('all')
        
# to plot v_los
#with sns.color_palette("RdPu",45):
#    for d in range(len(vf[0,0,0,:])):
#        plt.plot((-vf[0,0,:,d]).squeeze()*1e-5)
        
