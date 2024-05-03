import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
#import spectral as s
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


# DECLERATIONS
indx_case = 1
fibril_indxs = ['06_142705','06_165307','06_165721','07_130930'] # ['06_165307', '06_165307', '06_165307', '06_165307']
bg_indxs =  ['06_142809','06_165404','06_165752','07_131005'] # ['06_165404', '06_165404', '06_165404', '06_165404']
index = fibril_indxs[indx_case]
index_bg = bg_indxs[indx_case]
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

mode = ''
cyc = '1'
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

#settings board!
reg = 'n$_{T}$ = '+" ".join(str('%.1f'%(x)) for x in nodes_temp)+' \nn$_{v_{LOS}}$ = '+" ".join(str('%.1f'%(x)) for x in nodes_vlos)+' \nv$_{turb}$ = '+" ".join(str('%.1f'%(x)) for x in nodes_vturb)+'  \nreg$_{type}$ = '+" ".join(str(int(x)) for x in regularization_type)+'        reg = '+str('%.1f'%(regularize))+'       reg$_{weight}$ = '+" ".join(str('%.1f'%(x)) for x in regularization_weights)

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

# maps from the data
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
#file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
#file_ca8 = file_search(datadir,'crispex*'+pref[1]+'*.fcube')
file_cak = file_search(datadir,'crispex*'+pref[2]+'*.fcube')
#file_ha = file_search(datadir,'crispex*'+pref[3]+'*.fcube')
#cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])
#cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])
cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])
#cube_ha = lp_read(datadir+file_ha[0],datadir+file_ha[1])
c_map = cube_cak[fr,0,-1,:,:]
cak_int_un = (lp_read_scan(datadir+'cak_int_un.fcube'))[fr,:,:]

#calibration info
calib_cak = mf.readfits(datadir+'calib.3934.fits')
calib_ha = mf.readfits(datadir+'calib.6563.fits')

#inversion results
f_i = sp.profile(resdir+'fb/f'+mode+'_observed.nc')
f_o = sp.profile(resdir+'fb/f'+mode+'_synthetic_cycle'+cyc+'.nc')
f_m = sp.model(resdir+'fb/f'+mode+'_atmosout_cycle'+cyc+'.nc')
b_i = sp.profile(resdir+'bg/b'+mode+'_observed.nc')
b_o = sp.profile(resdir+'bg/b'+mode+'_synthetic_cycle'+cyc+'.nc')
b_m = sp.model(resdir+'bg/b'+mode+'_atmosout_cycle'+cyc+'.nc')

# Mapping the big array of results according to the individual fibril
ffile_fe = file_search(savedir, 'f_obs6302_*.fits')
ffile_ca8 = file_search(savedir, 'f_obs8542_*.fits')
ffile_ck = file_search(savedir, 'f_obs3950_*.fits')
bfile_fe = file_search(savedir, 'b_obs6302_*.fits')
bfile_ca8 = file_search(savedir, 'b_obs8542_*.fits')
bfile_ck = file_search(savedir, 'b_obs3950_*.fits')
#index = ffile_fe
finit = 0
binit = 0
I_dif = np.zeros(len(ffile_fe), dtype = float)
I_dif_h = np.zeros(len(ffile_fe), dtype = float)

for fn in range(len(ffile_fe)):
#for fn in range(10):
    fibril = mf.readfits(savedir + ffile_fe[fn])
    dark = mf.readfits(savedir + bfile_fe[fn])
    fnpix = fibril.shape[2]
    bnpix = dark.shape[2]

    #Fibril extract
    ##========
    # inversion
    norm = np.max(np.mean(f_o.dat[0,0,finit:finit+fnpix,:,ss] ,axis = 0))
    f_obs = 100*f_i.dat[0,0,finit:finit+fnpix,:,ss]/norm
    f_syn = 100*f_o.dat[0,0,finit:finit+fnpix,:,ss]/norm
    f_wav = f_i.wav
    f_n = fnpix
    f_dep = f_m.ltau[0,0,0,:]
    
    # slab
    fibdir = datadir+'fr'+str(fr)+'/'
    fib_file = (file_search(fibdir,'crispex*3950*.csav'))[2*fn]
    fib = restore(fibdir+fib_file)
    f_slab = fib.loop_slab[:,fr,:]*1000  #intensity values in the desired frame
    f_x = fib.x_coords
    f_x_pts = fib.x_loop_pts
    f_y = fib.y_coords
    f_l = fib.loop_size
    f_y_pts = fib.y_loop_pts
    
    # calculating the effective intensity
    # (i.e. the range of wl that contributes to chromospheric brightness)
    f_Itot = np.sum(f_slab[6:15,:])/f_slab[6:15,:].size

    #fibril in Ha
    fib_file_h = (file_search(fibdir,'crispex*6563*.csav'))[2*fn]
    fib_h = restore(fibdir+fib_file_h)
    f_slab_h = fib_h.loop_slab[:,fr,:]/calib_ha[0]  #intensity values in the desired frame
    
    # calculating the effective intensity
    # (i.e. the range of wl that contributes to chromospheric brightness)
    f_Itot_h = np.sum(f_slab_h[4:10,:])/f_slab_h[4:10,:].size

    #Background extract
    #===========
    #inversion
    b_obs = 100*b_i.dat[0,0,binit:binit+bnpix,:,ss]/norm
    b_syn = 100*b_o.dat[0,0,binit:binit+bnpix,:,ss]/norm
    b_wav = b_i.wav
    b_n = bnpix
    b_dep = b_m.ltau[0,0,0,:]
    
    # slab
    bg_file = file_search(fibdir,'crispex*3950*.csav')[2*fn+1]
    bg = restore(fibdir+bg_file)
    b_slab = bg.loop_slab[:,fr,:]*1000 #intensity values in the desired frame
    b_x = bg.x_coords
    b_x_pts = bg.x_loop_pts
    b_y_pts = bg.y_loop_pts
    b_y = bg.y_coords
    b_l = bg.loop_size
    
    # calculating the effective intensity
    # (i.e. the range of wl that contributes to chromospheric brightness)
    b_Itot = np.sum(b_slab[6:15,:])/b_slab[6:15,:].size

    # fibril/background valid length
    f_n_valid = np.min([f_n,b_n])
    n = 20
    tp = np.array([f_n_valid/n, f_n_valid/2, (n-1)*f_n_valid/n], dtype='uint')
    gg = len(tp)

    # fibril in Ha
    bg_file_h = (file_search(fibdir,'crispex*6563*.csav'))[2*fn+1]
    bg_h = restore(fibdir+bg_file_h)
    b_slab_h = bg_h.loop_slab[:,fr,:]/calib_ha[0]  #intensity values in the desired frame
    
    # calculating the effective intensity
    # (i.e. the range of wl that contributes to chromospheric brightness)
    b_Itot_h = np.sum(b_slab_h[4:10,:])/b_slab_h[4:10,:].size

    # FIBRIL and BG IDs:
    index = ((fib_file.split('_')[5]+'_'+fib_file.split('_')[6]).split('.')[0]).split('Aug')[1]
    index_bg = ((bg_file.split('_')[5]+'_'+bg_file.split('_')[6]).split('.')[0]).split('Aug')[1]
    for id in range(len(indxs)):
        if indxs[id] == index:
            #print 'f, b = ', [index,index_bg], 'id = ', fn
                    #print 'background: ',bg_file,'    ', bnpix, '=?', len(b_y_pts)

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
                            
            # wav = s.air2vac(i.wav)
            if(ss == 0):
                ymin = 0.2
                ymax = np.max(np.mean(f_obs, axis = 0))
            elif(ss == 1 or ss == 2):
                ymin = -0.055
                ymax = 0.005
            else:
                ymin = -0.03
                ymax = 0.03

            # amplifying the cak spectral profile
            amp = 1. #cak amplifier
            f_obs_amp = f_i.dat[0,0,finit:finit+fnpix,:,ss]
            f_obs_amp[:,0:42]  = f_i.dat[0,0,finit:finit+fnpix,0:42,ss]*amp
            f_syn_amp = f_o.dat[0,0,finit:finit+fnpix,:,ss]
            f_syn_amp[:,0:42]  = f_o.dat[0,0,finit:finit+fnpix,0:42,ss]*amp
            b_obs_amp = b_i.dat[0,0,binit:binit+bnpix,:,ss]
            b_obs_amp[:,0:42] = b_i.dat[0,0,binit:binit+bnpix,0:42,ss]*amp
            b_syn_amp = b_o.dat[0,0,binit:binit+bnpix,:,ss]
            b_syn_amp[:,0:42] = b_o.dat[0,0,binit:binit+bnpix,0:42,ss]*amp
            
            ##############
            #plot each fibril profile#
            ##############
            if plot_mode==1:
                plt.close("all")

                f = plt.figure(figsize=(8,7))
                #ax1 = plt.subplot2grid((10,6), (0,0), colspan=6, rowspan=3)
                ax11 = plt.subplot2grid((3,6), (1,0), colspan=2, rowspan=1)
                ax12 = plt.subplot2grid((3,6), (1,2), colspan=2, rowspan=1)
                ax13 = plt.subplot2grid((3,6), (1,4), colspan=2, rowspan=1)
                ax2 = plt.subplot2grid((3,6), (2,0), colspan=2, rowspan=1)
                ax4 = plt.subplot2grid((3,6), (2,2), colspan=2, rowspan=1)
                ax3 = plt.subplot2grid((3,6), (2,4), colspan=2, rowspan=1)
                ax5 = plt.subplot2grid((3,6),(0,0), colspan=3,rowspan=1)
                ax6 = plt.subplot2grid((3,6),(0,3), colspan=3,rowspan=1)
                
                line_i = '--'
                linewidth = 0.75
                markersize_i = 3
                ax11.plot(np.mean(b_obs,axis = 0)[0:42],'.',color = 'gray',markersize = 4.5, alpha = alp2)
                ax11.plot(np.mean(f_obs, axis = 0)[0:42],'.',color = 'red',markersize = 4.5, alpha = alp2)
                ax11.plot(np.mean(b_syn,axis = 0)[0:42], color='gray',linewidth = linewidth, alpha = alp2, linestyle = line_i)
                ax11.plot(np.mean(f_syn, axis = 0)[0:42], color='red',linewidth = linewidth, alpha = alp2, linestyle = line_i)

                ax12.plot(np.mean(b_obs,axis = 0)[43:129],'.',color = 'gray',markersize = 4.5)
                ax12.plot(np.mean(f_obs, axis = 0)[43:129],'.',color = 'red',markersize = 4.5)
                ax12.plot(np.mean(b_syn,axis = 0)[43:129], color='gray',linewidth = linewidth, alpha = alp2, linestyle = line_i)
                ax12.plot(np.mean(f_syn, axis = 0)[43:129], color='r',linewidth = linewidth, alpha = alp2, linestyle = line_i)

                ax13.plot(np.mean(b_obs,axis = 0)[130:],'.',color = 'gray',markersize = 4.5, label = 'b$\mathrm{_{obs}}$')
                ax13.plot(np.mean(f_obs, axis = 0)[130:],'.',color = 'red',markersize = 4.5, label = 'f$\mathrm{_{obs}}$')
                ax13.plot(np.mean(b_syn,axis = 0)[130:], color='gray',linewidth = linewidth, label = 'b$\mathrm{_{syn}}$', alpha = alp2, linestyle = line_i)
                ax13.plot(np.mean(f_syn, axis = 0)[130:], color='r',linewidth = linewidth, label = 'f$\mathrm{_{syn}}$', alpha = alp2, linestyle = line_i)

                
                for bb in range(0,b_n):
                    ax2.plot(b_m.temp[0,0,bb+binit,:].squeeze()*1e-3, b_dep, 'k-',color = 'gray',alpha = alp)
                    ax3.plot(b_m.vturb[0,0,bb+binit,:].squeeze()*1.e-5, b_dep, 'k-',color = 'gray',alpha = alp)
                    ax4.plot(b_m.vlos[0,0,bb+binit,:].squeeze()*1.e-5, b_dep, 'k-',color = 'gray',alpha = alp)
                for dd in range(45):
                    b_temp_m[dd] = np.mean(b_m.temp[0,0,binit:binit+bnpix,dd])
                    b_temp_med[dd] = np.median(b_m.temp[0,0,binit:binit+bnpix,dd])
                    b_los_m[dd] = np.mean(b_m.vlos[0,0,binit:binit+bnpix,dd])
                    b_los_med[dd] = np.median(b_m.vlos[0,0,binit:binit+bnpix,dd])
                    b_turb_m[dd] = np.mean(b_m.vturb[0,0,binit:binit+bnpix,dd])
                    b_turb_med[dd] = np.median(b_m.vturb[0,0,binit:binit+bnpix,dd])
                    
                for ff in range(0,f_n):
                    f_temp_m += f_m.temp[0,0,finit+ff,:]
                    ax2.plot(f_m.temp[0,0,finit+ff,:].squeeze()*1e-3,f_dep,  'k-',color = 'red',alpha = alp)
                    ax3.plot(f_m.vturb[0,0,finit+ff,:].squeeze()*1.e-5, f_dep, 'k-',color = 'red',alpha = alp)
                    ax4.plot( f_m.vlos[0,0,finit+ff,:].squeeze()*1.e-5, f_dep,'k-',color = 'red',alpha = alp)
                for dd in range(45):
                    f_temp_m[dd] = np.mean(f_m.temp[0,0,finit:finit+fnpix,dd])
                    f_temp_med[dd] = np.median(f_m.temp[0,0,finit:finit+fnpix,dd])
                    f_los_m[dd] = np.mean(f_m.vlos[0,0,finit:finit+fnpix,dd])
                    f_los_med[dd] = np.median(f_m.vlos[0,0,finit:finit+fnpix,dd])
                    f_turb_m[dd] = np.mean(f_m.vturb[0,0,finit:finit+fnpix,dd])
                    f_turb_med[dd] = np.median(f_m.vturb[0,0,finit:finit+fnpix,dd])
                    
                lw = 1.5 #plot linewidth
                ax2.plot((b_temp_m).squeeze()*1e-3, b_dep, 'k-',color = 'gray',linewidth = lw, alpha = alp2, label = 'b')
                ax2.plot((f_temp_m).squeeze()*1e-3, f_dep, 'k-',color = 'r',linewidth = lw, alpha = alp2,  label = 'f')
                for n in range(len(nodes_temp)):
                    ax2.axhline(nodes_temp[n],linestyle = '--',color = 'k',linewidth = 0.1, alpha = 0.6)
                    ax3.plot((b_turb_m).squeeze()*1e-5, b_dep, 'k-',color = 'gray',linewidth = lw, alpha = alp2)
                    ax3.plot((f_turb_m).squeeze()*1e-5, f_dep, 'k-',color = 'r',linewidth = lw, alpha = alp2)
                for n in range(len(nodes_vturb)):
                    ax3.axhline(nodes_vturb[n],linestyle = '--',color = 'k',linewidth = 0.1, alpha = 0.6)
                    ax4.plot((b_los_m).squeeze()*1e-5, b_dep, 'k-',color = 'gray',linewidth = lw, alpha = alp2)
                    ax4.plot((f_los_m).squeeze()*1e-5, f_dep, 'k-',color = 'r',linewidth = lw, alpha = alp2)
                for n in range(len(nodes_vlos)):
                    ax4.axhline(nodes_vlos[n],linestyle = '--',color = 'k',linewidth = 0.1, alpha = 0.6)
                
                    ax11.set_ylabel(r'$I$ / $I_{\mathrm{cont}}$ [%]', fontdict = font)
                    ax11.set_xlabel(r'$\lambda$ [$\mathrm{\AA}$]', fontdict = font)
                    cak_tick_pos = [10,21,32]
                    ax11.set_xticks(cak_tick_pos)
                    ax11.set_xticklabels(np.round(f_i.wav[cak_tick_pos], decimals = 1), fontdict = font)
                    ax11.set_ylim(7,18)
                    ax11.set_xlim(0,41)
                    ax11.tick_params(axis='both', which='major', labelsize=8)

                    ax12.set_xlabel(r'$\lambda$ [$\mathrm{\AA}$]', fontdict = font)
                    ca8_tick_pos = [86-43-22, 86-43, 86-43+22]
                    ax12.set_xticks(ca8_tick_pos)
                    ax12.set_xticklabels(np.round(f_i.wav[np.array(ca8_tick_pos)+43], decimals = 1), fontdict = font)
                    ax12.set_ylim(20,72)
                    ax12.set_xlim(0,128-43)
                    ax12.tick_params(axis='both', which='major', labelsize=8)

                    ax13.set_xlabel(r'$\lambda$ [$\mathrm{\AA}$]', fontdict = font)
                    fe_tick_pos = [153 -129,182-129,211-129]
                    ax13.set_xticks(fe_tick_pos)
                    ax13.set_xticklabels(np.round(f_i.wav[np.array(fe_tick_pos)+129], decimals = 1), fontdict = font)
                    ax13.set_ylim(40,100)
                    ax13.set_xlim(1,b_wav.size-1-129)
                    f.legend(ncol = 6, fontsize = 7, bbox_to_anchor=(0.25, 0.185, 0.5, 0.5), framealpha = 0.2)
                    ax13.tick_params(axis='both', which='major', labelsize=8)

    
                ax2.set_xlim(3.75,8.5)
                ax2.set_ylim(0.1,-6.3)
                ax2.set_ylabel(r'log($\tau_{500}$)', fontdict = font)
                ax2.set_xlabel(r'$T$ [kK]', fontdict = font)
                ax2.tick_params(axis='both', which='major', labelsize=8)

                ax3.set_xlim(-0.1,8.25)
                ax3.set_ylim(0.1,-6.3)
                #ax3.set_ylabel(r'log($\mathrm{\tau}$)', fontdict = font)
                ax3.set_xlabel(r'$v_{\mathrm{turb}}$ [km s$^{-1}$]', fontdict = font)
                ax3.tick_params(axis='both', which='major', labelsize=8)
                                
                ax4.set_xlim(-6,7.5)
                ax4.set_ylim(0.1,-6.3)
                #ax4.set_ylabel(r'log($\mathrm{\tau}$)', fontdict = font)
                ax4.set_xlabel(r'$v\mathrm{_{LOS}}$ [km s$^{-1}$]', fontdict = font)
                ax4.tick_params(axis='both', which='major', labelsize=8)

                #cropping maps
                edge = 40
                xmin = int(np.min([f_x_pts[0],f_x_pts[-1]])-edge)
                xmax = int(np.max([f_x_pts[0],f_x_pts[-1]])+edge)
                ymin = int(np.min([f_y_pts[0],f_y_pts[-1]])-edge)
                ymax = int(np.max([f_y_pts[0],f_y_pts[-1]])+edge)
                xlabel = 'x [arcsec]'
                ylabel = 'y [arcsec]'
                cmap = 'gray'
                title = r'Continuum 4000 $\mathrm{\AA}$'
                titlek = r'Ca II K $\mathrm{\lambda}$-integrated'
                ytick_pos = np.arange(0,(np.round((ymax-ymin)/40)+1)*40,40)
                ytick_lab = ytick_pos*res
                xtick_pos = np.arange(0,(np.round((xmax-xmin)/40)+1)*40,40)
                xtick_lab = xtick_pos*res
                
                #cont map
                c_map_cropped = c_map[ymin:ymax,xmin:xmax]
                ax5.imshow(c_map_cropped, cmap = cmap)
                ax5.set_xlim(0, xmax-xmin)
                ax5.set_ylim(0, ymax-ymin)
                ax5.set_xlabel(xlabel, fontdict = font)
                ax5.set_ylabel(ylabel, fontdict = font)
                ax5.xaxis.set_minor_locator(AutoMinorLocator(5))
                ax5.yaxis.set_minor_locator(AutoMinorLocator(5))
                ax5.tick_params(which='minor', length=2)
                ax5.set_yticks(ytick_pos)
                ax5.set_yticklabels(ytick_lab, fontdict=font)
                ax5.set_xticks(xtick_pos)
                ax5.set_xticklabels(xtick_lab, fontdict=font)

                #ax5.set_title(title)
                ax5.text(1,0.5,r'Continuum 4000 $\mathrm{\AA}$',color = 'white',horizontalalignment='left', verticalalignment='bottom', fontdict = font) #, bbox={'facecolor': 'black', 'pad': 10,'alpha':0.5}
                
                # overplotting paths on the maps
                ax5.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'gray',alpha = 0.85)
                # overplotting the fibril
                ax5.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'red',alpha = 0.85)
                
                # cak map
                cak_map_cropped = cak_int_un[ymin:ymax,xmin:xmax]
                ax6.imshow(cak_map_cropped, cmap = cmap)
                ax6.set_xlim(0, xmax-xmin)
                ax6.set_ylim(0, ymax-ymin)
                ax6.set_xlabel(xlabel, fontdict = font)
                ax6.xaxis.set_minor_locator(AutoMinorLocator(5))
                ax6.yaxis.set_minor_locator(AutoMinorLocator(5))
                ax6.tick_params(which='minor', length=2)
                ax6.set_yticks(ytick_pos)
                ax6.set_yticklabels(ytick_lab, fontdict=font)
                ax6.set_xticks(xtick_pos)
                ax6.set_xticklabels(xtick_lab, fontdict=font)
                ax6.set_ylabel('')
                ax6.text(1,0.5,r'Ca II K $\mathrm{\lambda}$-integrated',color = 'white',horizontalalignment='left', verticalalignment='bottom', fontdict = font) #, bbox={'facecolor': 'black', 'pad': 10,'alpha':0.5}
                # overplotting the fibril on the maps
                ax6.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'gray',alpha = 0.4,linewidth = 4)
                ax6.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'red',alpha = 0.3,linewidth = 4)

                for ii in range(gg):
                    ax5.plot(b_x_pts[tp[ii]]-xmin,b_y_pts[tp[ii]]-ymin, marker = markers[ii], color = 'gray',markersize = 5, alpha = 0.5)
                    ax5.plot(f_x_pts[tp[ii]]-xmin,f_y_pts[tp[ii]]-ymin, marker = markers[ii], color = 'red',markersize = 5, alpha = 0.5)
                    ax6.plot(b_x_pts[tp[ii]]-xmin,b_y_pts[tp[ii]]-ymin, marker = markers[ii], color = 'gray',markersize = 5,markerfacecolor = 'none')
                    ax6.plot(f_x_pts[tp[ii]]-xmin,f_y_pts[tp[ii]]-ymin, marker = markers[ii], color = 'red',markersize = 5,markerfacecolor = 'none')

                f.set_tight_layout(True)
                plt.ion()
                plt.show()

                outname = outdir+'inv_pro_'+index+'_'+mode+'_cycle'+cyc+'.pdf'
                f.savefig(outname)
                print 'file saved to:'+outname

    if fn ==95:
        break

    finit = finit + fnpix
    binit = binit + bnpix
