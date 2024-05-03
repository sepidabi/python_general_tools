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

# fibril
index = '06_165307'
index_bg = '06_165404'
tp = [27, 123, 210] # pixels of interest along the fibril
gg = len(tp) # number of pixels of interest

# plotting font, etc.
plot_mode = 3    # set to 1 for fibril prof
                             # 2 for RF
                             # 3 for along the fibril
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }
markers = ['s','o','^'] # for the desired pixels
alp = 0.01 #transparency of the plotted lines

# inversion info
mode = 'RF' # set to 'RF' for the response function (grid) 
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

datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
outdir = '/home/seki2695/OUTPUT/inv/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
fn = 28 #frame no.
ss = 0 # stokes param index
#xx = 100  # set to zero to choose the first pixel
#yy = 0     #stokes I values in the dat array

b_i = sp.profile('/home/seki2695/INV/stic/bg/b'+mode+'_observed_'+index+'.nc')
b_o = sp.profile('/home/seki2695/INV/stic/bg/b'+mode+'_synthetic_cycle'+cyc+'_'+index+'.nc')
b_m = sp.model('/home/seki2695/INV/stic/bg/b'+mode+'_atmosout_cycle'+cyc+'_'+index+'.nc')
b_mRF = sp.model('/home/seki2695/INV/stic/bg/b'+mode+'_atmosout_cycle'+cyc+'_'+index+'.nc')
b_rf = sp.profile('/home/seki2695/INV/stic/bg/b'+mode+'_synthetic_cycle2_'+index+'.nc')
#b_i = sp.profile('/home/seki2695/INV/stic/bg/b'+mode+'_observed_'+index+'.nc')
#b_o = sp.profile('/home/seki2695/INV/stic/bg/b'+mode+'_synthetic_cycle'+cyc+'_'+index+'.nc')
#b_m = sp.model('/home/seki2695/INV/stic/bg/b'+mode+'_atmosout_cycle'+cyc+'_'+index+'.nc')


#Fibril extract
##========
amp = 1. #cak amplifier
#inversion
b_obs = b_i.dat[0,0,:,:,ss]
b_obs[:,0:42]  = b_i.dat[0,0,:,0:42,ss]*amp
b_syn = b_o.dat[0,0,:,:,ss]
b_syn[:,0:42]  = b_o.dat[0,0,:,0:42,ss]*amp
b_wav = b_i.wav
b_n = b_i.dat.shape[2]
b_dep = b_m.ltau[0,0,0,:]
b_depRF = b_mRF.ltau[0,0,0,:]
#slab
fibdir = datadir+'fr'+str(fn)+'/'
fib_file = file_search(fibdir,'crispex*3950*'+index+'*.csav')
fib = restore(fibdir+fib_file[0])
b_slab = fib.loop_slab  #intensity values in the desired frame
b_x = fib.x_coords
b_x_pts = fib.x_loop_pts
b_y = fib.y_coords
b_l = fib.loop_size
b_y_pts = fib.y_loop_pts


#containers:
b_temp_m = np.zeros(45)
b_temp_med = np.zeros(45)
b_los_m = np.zeros(45)
b_los_med = np.zeros(45)
b_turb_m = np.zeros(45)
b_turb_med = np.zeros(45)

#chi2
pxnl = np.min([b_n])
chi2_f = np.zeros(pxnl)

for px in range(0,pxnl):
    chi2_f[px] = np.sum(((b_o.dat[0,0,px,:,ss] - b_i.dat[0,0,px,:,ss])/b_o.weights[:,ss])**2)/len(b_wav)
chi2_fm = np.mean(chi2_f)

#maps
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

#fibril
fibdir = datadir+'fr'+str(fn)+'/'
fib_file = file_search(fibdir,'crispex*3950*'+index+'*.csav')
fib = restore(fibdir+fib_file[0])
b_slab = fib.loop_slab  #intensity values in the desired frame
b_x = fib.x_coords
b_x_pts = fib.x_loop_pts
b_y = fib.y_coords
b_l = fib.loop_size
b_y_pts = fib.y_loop_pts

#wav = s.air2vac(i.wav)
if(ss == 0):
    ymin = 0.05
    ymax = np.max(np.mean(b_obs, axis = 0))
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
    b_obs_amp = b_i.dat[0,0,:,:,ss]
    b_obs_amp[:,0:42]  = b_i.dat[0,0,:,0:42,ss]*amp
    b_syn_amp = b_o.dat[0,0,:,:,ss]
    b_syn_amp[:,0:42]  = b_o.dat[0,0,:,0:42,ss]*amp

    #ax1.plot(np.mean(b_obs,axis = 0),'s',color = 'blue',markersize = 3)
    ax1.plot(np.mean(b_obs, axis = 0),'.',color = 'red',markersize = 4.5)
    #ax1.plot(np.mean(b_syn,axis = 0), color='b',linewidth = 1.25)
    ax1.plot(np.mean(b_syn, axis = 0), color='r',linewidth = 1.25)
    
    for ff in range(0,b_n):
        b_temp_m += b_m.temp[0,0,ff,:]
        ax2.plot(b_dep, b_m.temp[0,0,ff,:].squeeze()*1e-3, 'k-',color = 'red',alpha = alp)
        ax3.plot(b_dep, b_m.vturb[0,0,ff,:].squeeze()*1.e-5, 'k-',color = 'red',alpha = alp)
        ax4.plot(b_dep, b_m.vlos[0,0,ff,:].squeeze()*1.e-5, 'k-',color = 'red',alpha = alp)
    for dd in range(45):
        b_temp_m[dd] = np.mean(b_m.temp[0,0,:,dd])
        b_temp_med[dd] = np.median(b_m.temp[0,0,:,dd])
        b_los_m[dd] = np.mean(b_m.vlos[0,0,:,dd])
        b_los_med[dd] = np.median(b_m.vlos[0,0,:,dd])
        b_turb_m[dd] = np.mean(b_m.vturb[0,0,:,dd])
        b_turb_med[dd] = np.median(b_m.vturb[0,0,:,dd])
    
    lw = 1.1 #plot linewidth
    #ax2.plot(b_dep, (b_temp_m).squeeze()*1e-3, 'k-',color = 'b',linewidth = lw)
    #ax2.plot(b_dep, (b_temp_med).squeeze()*1e-3, 'k-',color = 'b',linestyle = '--',linewidth = lw)
    ax2.plot(b_dep, (b_temp_m).squeeze()*1e-3, 'k-',color = 'r',linewidth = lw)
    #ax2.plot(b_dep, (b_temp_med).squeeze()*1e-3, 'k-',color = 'r',linestyle = '--',linewidth = lw)
    for n in range(len(nodes_temp)):
        ax2.axvline(nodes_temp[n],linestyle = '--',color = 'k',linewidth = 0.2)
        #ax3.plot(b_dep, (b_turb_m).squeeze()*1e-5, 'k-',color = 'b',linewidth = lw)
        #ax3.plot(b_dep, (b_turb_med).squeeze()*1e-5, 'k-',color = 'b',linestyle = '--',linewidth = lw)
        ax3.plot(b_dep, (b_turb_m).squeeze()*1e-5, 'k-',color = 'r',linewidth = lw)
        #ax3.plot(b_dep, (b_turb_med).squeeze()*1e-5, 'k-',color = 'r',linestyle = '--',linewidth = lw)
    for n in range(len(nodes_vturb)):
        ax3.axvline(nodes_vturb[n],linestyle = '--',color = 'k',linewidth = 0.2)
        #ax4.plot(b_dep, (b_los_m).squeeze()*1e-5, 'k-',color = 'b',linewidth = lw)
        #ax4.plot(b_dep, (b_los_med).squeeze()*1e-5, 'k-',color = 'b',linestyle = '--',linewidth = lw)
        ax4.plot(b_dep, (b_los_m).squeeze()*1e-5, 'k-',color = 'r',linewidth = lw)
        #ax4.plot(b_dep, (b_los_med).squeeze()*1e-5, 'k-',color = 'r',linestyle = '--',linewidth = lw)
    for n in range(len(nodes_vlos)):
        ax4.axvline(nodes_vlos[n],linestyle = '--',color = 'k',linewidth = 0.2)

    ax1.set_ylabel(r'I$_{mean}$ [erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$]')
    ax1.set_xlabel(r'Wavelength [$\mathrm{\AA}$]')
    ax1.set_xticks(([21,86,153,211]))
    ax1.set_xticklabels((['3934','8542','6301','6302']))
    ax1.set_ylim(ymin, ymax+0.05)
    #ax1.set_xlim(0,b_wav.size-1)
    ax1.legend(['Background','Fibril'])
    
    ax2.set_ylim(3.75,11)
    ax2.set_xlim(-6.5,1)
    ax2.set_xlabel(r'log($\mathrm{\tau}$)')
    ax2.set_ylabel(r'T [kK]')

    ax3.set_ylim(-0.1,10)
    ax3.set_xlim(-7.5,1)
    ax3.set_xlabel(r'log($\mathrm{\tau}$)')
    ax3.set_ylabel(r'v$_{\rm turb}$ [km/s]')

    ax4.set_ylim(-7.5,7.5)
    ax4.set_xlim(-7.5,1)
    ax4.set_xlabel(r'log($\mathrm{\tau}$)')
    ax4.set_ylabel(r'v$\mathrm{_{LOS}}$ [km/s]')

    #inv = restore('inv_in'+index+'.sav')
    #xmin = 0
    #xmax = c_map.shape[1]
    #ymin = 0
    #ymax = c_map.shape[0]
    xmin = int(np.min([b_x_pts[0],b_x_pts[-1]])-20)
    xmax = int(np.max([b_x_pts[0],b_x_pts[-1]])+20)
    ymin = int(np.min([b_y_pts[0],b_y_pts[-1]])-20)
    ymax = int(np.max([b_y_pts[0],b_y_pts[-1]])+20)
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

    #overplotting the background
    #ax5.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue')
    #ax5.plot(inv.b_x_pts[inv.pxn],inv.b_y_pts[inv.pxn], '+', color = 'blue',markersize = 15)
    #overplotting the fibril
    ax5.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'red')
    #ax5.plot(inv.b_x_pts[inv.pxn],inv.b_y_pts[inv.pxn], '+', color = 'red',markersize = 15)
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
    #ax6.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue',alpha = 0.3,linewidth = 4)
    ax6.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'red',alpha = 0.3,linewidth = 4)
    f.set_tight_layout(True)
    plt.ion()
    plt.show()

    f.savefig(outdir+'inv_pro_'+index+'_'+mode+'_cycle'+cyc+'.pdf')
    print 'file saved to:'+outdir+'inv_pro_'+index+'_'+mode+'_cycle'+cyc+'.pdf'

    #plt.close('all')
    
# Plotting RF #
##########
if plot_mode == 2:

    tp = [0,50,150] # pixels of interest index
    gg = len(tp) # number of pixels of interest

    plt.close('all')
    f = plt.figure(figsize = [3*gg-1.,3*gg])

    gs1 = gridspec.GridSpec(3, 1) # grid scale of the 1st col.
    gs1.update(left=0.03,
               right=0.33,
               #wspace=0.05,
               bottom=0.05,
               top=0.94,
               hspace = 0.5)
    gs = gridspec.GridSpec(3, 2) # grid scale of the 2nd & 3rd col.
    gs.update(left=0.41,
               right=1.,
               wspace = 0.05,
               bottom=0.05,
               top=0.94,
               hspace=0.5)

    for pp in range(gg):

        # Temperature RF
        ax1 = host_subplot(gs1[pp, 0], adjustable='box', aspect='equal')
        rb_t = b_rf.rf[0,0,:,0,:,0:143,0]
        norm_t = np.max(np.abs(rb_t[tp,:,:]))
        temperature = ax1.imshow(im.histo_opt(rb_t[tp[pp],:,:]),cmap = 'bone', vmin = 0, vmax = norm_t)
        #plt.subplot(2,1,1)
        #plt.colorbar()
        ax1.plot(np.argmax(rb_t[tp[pp],:,:], axis = 0), color = 'white',linewidth = 0.75)
        ax1.set_xlabel(r'Wavelength [$\mathrm{\AA}$]', fontdict=font)
        ax1.set_ylabel(r'log($\mathrm{\tau}$)', fontdict=font)
        xtick_pos = [14,72,130]
        xtick_lab = np.round(b_wav[xtick_pos], decimals = 1)
        ax1.set_xticks(xtick_pos)
        ax1.set_xticklabels(xtick_lab, fontdict=font)
        ytick_pos = [13, 30, 47, 64, 81, 98, 115, 132]#[8,18,28,38,48,58,68,78,88]
        ytick_lab = np.round(b_depRF[ytick_pos])
        ax1.set_yticks(ytick_pos)
        ax1.set_yticklabels(ytick_lab.astype(int), fontdict=font)
        
        ax12 = ax1.twinx()
        ax12.set_ylim(0.01, 0.4)
        ax12.set_ylabel(r'I [erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$]', fontdict=font)
        ax12.plot(b_obs[tp[pp],0:142],'.',color = 'b',marker = markers[pp], markersize = 3.)#, label = r'I$_{\rm obs}$')
        ax12.plot(b_syn[tp[pp],0:142], color='b',linewidth = 1., label = r'I$_{\rmsyn}$')
        
        ax13 = ax1.twiny()
        b_tempRF =b_mRF.temp[0,0,tp[pp],:]
        ax13y = np.linspace(0,len(b_depRF)-1,147)
        ax13.set_xlim(3.75,9)
        ax13.set_xlabel(r'T [kK]', fontdict=font)
        ax13.plot(b_tempRF.squeeze()*1e-3,ax13y, 'k-',color = 'steelblue',linewidth = 1., label = 'T')
        
        #Required to remove some white border               
        ax1.autoscale(axis = 'both', enable = True)

        # Vlos RF
        ax2 = host_subplot(gs[pp,0], adjustable='box', aspect='equal')
        rb_vlos = b_rf.rf[0,0,:,1,:,0:143,0]
        norm_vlos = np.max(np.abs(rb_vlos))
        vlos = ax2.imshow(im.histo_opt(rb_vlos[tp[pp],:,:]),cmap = 'bwr',aspect = 'equal',vmin = -norm_vlos, vmax = norm_vlos)
        ax2.set_xlabel(r'Wavelength [$\mathrm{\AA}$]', fontdict=font)
        ax2.set_ylabel(r'log($\mathrm{\tau}$)', fontdict=font)
        xtick_pos = [14,72,130]
        xtick_lab = np.round(b_wav[xtick_pos], decimals = 1)
        ax2.set_xticks(xtick_pos)
        ax2.set_xticklabels(xtick_lab, fontdict=font)
        ytick_pos = [13, 30, 47, 64, 81, 98, 115, 132]#[8,18,28,38,48,58,68,78,88]
        ytick_lab = np.round(b_depRF[ytick_pos])
        ax2.set_yticks(ytick_pos)
        ax2.set_yticklabels(ytick_lab.astype(int), fontdict=font)
        
        ax23 = ax2.twiny()
        b_vlosRF =b_mRF.vlos[0,0,tp[pp],:]
        ax23y = np.linspace(0,len(b_depRF)-1,147)
        ax23.set_xlim(-7.5,7.5)
        ax23.set_xlabel(r'v$\mathrm{_{LOS}}$ [km/s]', fontdict=font)
        ax23.plot(b_vlosRF.squeeze()*1e-5,ax23y, 'k-',color = 'black',linewidth = 1., label = r'v$\rm _{LOS}$')
        
        #Required to remove some white border               
        ax2.autoscale(axis = 'both', enable = True)

        
        # Vmt RF
        ax3 = host_subplot(gs[pp,1], adjustable='box', aspect='equal')
        rb_vmt = b_rf.rf[0,0,:,2,:,0:143,0]
        norm_vmt = np.max(np.abs(rb_vmt))
        vmt = ax3.imshow(im.histo_opt(rb_vmt[tp[pp],:,:]),cmap = 'bwr',aspect = 'equal', vmin = -norm_vmt, vmax = norm_vmt)
        ax3.set_xlabel(r'Wavelength [$\mathrm{\AA}$]', fontdict=font)
        ax3.set_ylabel(r'log($\mathrm{\tau}$)', fontdict=font)
        xtick_pos = [14,72,130]
        xtick_lab = np.round(b_wav[xtick_pos], decimals = 1)
        ax3.set_xticks(xtick_pos)
        ax3.set_xticklabels(xtick_lab, fontdict=font)
        ytick_pos = [13, 30, 47, 64, 81, 98, 115, 132]#[8,18,28,38,48,58,68,78,88]
        ytick_lab = np.round(b_depRF[ytick_pos])
        ax3.set_yticks(ytick_pos)
        ax3.set_yticklabels(ytick_lab.astype(int), fontdict=font)
               
        ax33 = ax3.twiny()
        b_vmtRF =b_mRF.vturb[0,0,tp[pp],:]
        ax33y = np.linspace(0,len(b_depRF)-1,147)
        ax33.set_xlim(-1,8)
        ax33.set_xlabel(r'v$_{\rm turb}$ [km/s]', fontdict=font)
        ax33.plot(b_vmtRF.squeeze()*1e-5,ax33y, 'k-',color = 'black',linewidth = 1.,linestyle='--', label = r'v$_{\rm turb}$')
        
        #Required to remove some white border               
        ax3.autoscale(axis = 'both', enable = True)

        if pp==0:
            ax1.legend(loc = 'upper left')
            ax2.legend()
            ax3.legend()
            
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

    #f.tight_layout()
    plt.show()
    f.savefig(outdir+'inv_RF_'+index+'_'+mode+'_cycle'+cyc+'bg.eps',quality = 100)
    print 'file saved to:'+outdir+'inv_RF_'+index+'_'+mode+'_cycle'+cyc+'bg.eps'

################
# Variations along the fibril
################
if plot_mode ==3:

    plt.close('all')

    f = plt.figure(figsize = [5.5,6])
    gs = gridspec.GridSpec(3, 1) # grid scale of the 1st col.
    aspect = 0.775 # panels aspect ratio (y/x)
    gs.update(left=0.075,
               right=0.88,
               #wspace=0.05,
               bottom=0.08,
               top=0.99,
               hspace = 0.0)

    ax1 = plt.subplot(gs[0,0],adjustable = 'box',aspect = 'equal')
    ax2 = plt.subplot(gs[1,0],adjustable = 'box',aspect = 'equal')
    ax3 = plt.subplot(gs[2,0],adjustable = 'box',aspect = 'equal')

    # clipping the depth between log = [0,-7]
    dep_init = 13 # logtau = -7
    dep_fin = 132 # logtau = 0
    b_depRF = b_depRF[dep_init: dep_fin+1]
    depth_grid_sp = 17

    # normal depth grids
    #ytick_pos = [4, 9, 14, 19, 24, 29, 34, 39, 44]#[8,18,28,38,48,58,68,78,88]
    #ytick_lab = np.round(b_dep[ytick_pos])
    # RF depth grids
    ytick_pos = np.linspace(1,7,7,dtype='uint')*depth_grid_sp
    ytick_lab = np.round(b_depRF[ytick_pos])


    # temperature
    temperature = np.transpose(b_m.temp[0,0,:,dep_init:dep_fin+1].squeeze()*1e-3)
    #norm_t = np.max(np.abs(temperature))
    panel1 = ax1.imshow(temperature, cmap = 'bone', vmin = 3.75, vmax = 11,aspect = aspect)
    #ax1.set_xlabel(r'Wavelength [$\mathrm{\AA}$]', fontdict=font)
    ax1.set_ylabel(r'log($\mathrm{\tau}$)', fontdict=font)
    ax1.set_xticklabels([])
    ax1.set_yticks(ytick_pos)
    ax1.set_yticklabels(ytick_lab.astype(int), fontdict=font)
    for tt in range(gg):
        ax1.axvline(tp[tt], color = 'white', linestyle = '--')
        ax2.axvline(tp[tt], linestyle = '--', color = 'black')
        ax3.axvline(tp[tt], linestyle = '--', color = 'black')

    # temp colorbar
    axins = inset_axes(ax1, #axis
                       width="2%",  # width = 10% of parent_bbox width
                       height="90%",  # height : 50%
                       loc='center left',
                       bbox_to_anchor=(1.01, 0., 1, 1),
                       bbox_transform=ax1.transAxes,
                       borderpad=0,
    )
    ct_ticks = np.linspace(4,11,8,dtype='uint')
    ct = plt.colorbar(panel1, cax=axins, orientation="vertical",
                      ticks = ct_ticks,
    )
    ct.set_label(r'T (kK)', fontdict = font)
    ct.ax.set_yticklabels(['4', '5', '6', '7', '8', '9', '10', '$>$11'], fontdict=font)

    # vlos
    vlos = np.transpose(b_m.vlos[0,0,:,dep_init:dep_fin+1].squeeze()*1e-5)
    norm_v = np.max(np.abs(vlos))
    panel2 = ax2.imshow(vlos, cmap = 'bwr', vmin = -norm_v, vmax = norm_v,aspect = aspect)
    ax2.set_ylabel(r'log($\mathrm{\tau}$)', fontdict=font)
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
    vmt = np.transpose(b_m.vturb[0,0,:,dep_init:dep_fin+1].squeeze()*1e-5)
    norm_vmt = np.max(np.abs(vmt))
    panel3 = ax3.imshow(vmt, cmap = 'Reds', vmin = 0, vmax = norm_vmt,aspect = aspect)
    ax3.set_xlabel(r'pixel along the fibril', fontdict=font)
    ax3.set_ylabel(r'log($\mathrm{\tau}$)', fontdict=font)
    ax3.set_yticks(ytick_pos)
    ax3.set_yticklabels(ytick_lab.astype(int), fontdict=font)

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
    cvmt.set_label(r'v$_{\rm turb}$ [km/s]', fontdict = font)

    plt.show()

    f.savefig(outdir+'inv_pro_'+index+'_'+mode+'_cycle'+cyc+'lengthbg.eps', quality = 100)
    print 'file saved to:'+outdir+'inv_pro_'+index+'_'+mode+'_cycle'+cyc+'lengthbg.eps'


    #plt.close('all')

# to plot v_los
#with sns.color_palette("RdPu",45):
#    for d in range(len(vf[0,0,0,:])):
#        plt.plot((-vf[0,0,:,d]).squeeze()*1e-5)
