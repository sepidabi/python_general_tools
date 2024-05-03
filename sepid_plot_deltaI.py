import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
from matplotlib.gridspec import GridSpec

fontsize = 10.
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': fontsize,
}


# DECLERATIONS
indxs = ['06_165307', '06_142705','06_165721','07_130930']

plot_mode = 3 # 1 for each fibril profile
                          # 2 for (manual) dI histogram
                          # 3 for seaborn d(I) scatter plot
                          # 4 for seaborn d(T) scatter plot

mode = '_fullstokes'
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

alp = 0.01 #transparency of the plotted lines
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
savedir = datadir+'OUTPUT/'
resdir = '/home/seki2695/INV/stic/'
outdir = '/home/seki2695/OUTPUT/inv/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
fr = 28 #frame no.
fibdir = datadir+'fr'+str(fr)+'/'
ss = 0 # stokes param index
#xx = 100  # set to zero to choose the first pixel
#yy = 0     #stokes I values in the dat array

#maps from the data
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

# calibration info
calib_cak = mf.readfits(datadir+'calib.3934.fits')
calib_ha = mf.readfits(datadir+'calib.6563.fits')
calib_ca8 = mf.readfits(datadir+'calib.8542.fits')

# inversion results
f_i = sp.profile(resdir+'fb_back/f'+mode+'_observed.nc')
f_o = sp.profile(resdir+'fb_back/f'+mode+'_synthetic_cycle'+cyc+'.nc')
f_m = sp.model(resdir+'fb_back/f'+mode+'_atmosout_cycle'+cyc+'.nc')
b_i = sp.profile(resdir+'bg_back/b'+mode+'_observed.nc')
b_o = sp.profile(resdir+'bg_back/b'+mode+'_synthetic_cycle'+cyc+'.nc')
b_m = sp.model(resdir+'bg_back/b'+mode+'_atmosout_cycle'+cyc+'.nc')
# number of paired f-b pixels
f_b_pair_n = len(b_m.temp[0,0,:,0])

#Mapping the big array of results according to the individual fibril
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
I_dif_ca8 = np.zeros(len(ffile_fe), dtype = float)
I_dif_each = np.zeros(f_b_pair_n, dtype = float)
t_dif_mean = np.zeros(len(ffile_fe), dtype = float)
t_dif_each = np.zeros(f_b_pair_n, dtype = float)
approx_cak_dep_n = 28 # corresponding to log(tau) = -4.59 (RFmax)

for fn in range(len(file_search(fibdir,'crispex*3950*_2018*.csav'))/2):
#for fn in range(10):
    fibril = mf.readfits(savedir + ffile_fe[fn])
    dark = mf.readfits(savedir + bfile_fe[fn])
    fnpix = fibril.shape[2]
    bnpix = dark.shape[2]

    #Fibril extract
    ##========
    amp = 1. #cak amplifier
    #inversion
    f_obs = f_i.dat[0,0,finit:finit+fnpix,:,ss]
    f_obs[:,0:42]  = f_i.dat[0,0,finit:finit+fnpix,0:42,ss]#*amp
    f_syn = f_o.dat[0,0,finit:finit+fnpix,:,ss]
    f_syn[:,0:42]  = f_o.dat[0,0,finit:finit+fnpix,0:42,ss]#*amp
    f_wav = f_i.wav
    f_n = fnpix
    f_dep = f_m.ltau[0,0,0,:]
    ndep = len(f_dep)
    #slab
    fib_file = (file_search(fibdir,'crispex*3950*_2018*.csav'))[2*fn]
    fib = restore(fibdir+fib_file)
    f_slab = fib.loop_slab[:,fr,:]*1000  #intensity values in the desired frame
    f_x = fib.x_coords
    f_x_pts = fib.x_loop_pts
    f_y = fib.y_coords
    f_l = fib.loop_size
    f_y_pts = fib.y_loop_pts
    #calculating the effective intensity
    #(i.e. the range of wl that contributes to chromospheric brightness)
    f_Itot = np.sum(f_slab[6:15,:])/f_slab[6:15,:].size

    #fibril in Ha
    fib_file_h = (file_search(fibdir,'crispex*6563*_2018*.csav'))[2*fn]
    fib_h = restore(fibdir+fib_file_h)
    f_slab_h = fib_h.loop_slab[:,fr,:]/calib_ha[0]  #intensity values in the desired frame
    #calculating the effective intensity
    #(i.e. the range of wl that contributes to chromospheric brightness)
    f_Itot_h = np.sum(f_slab_h[4:10,:])/f_slab_h[4:10,:].size
    
    #fibril in Ca 8542
    fib_file_ca8 = (file_search(fibdir,'crispex*8542*_2018*.csav'))[2*fn]
    fib_ca8 = restore(fibdir+fib_file_ca8)
    f_slab_ca8 = fib_ca8.loop_slab[0,:,fr,:]/calib_ca8[0]  #intensity values in the desired frame
    #calculating the effective intensity
    #(i.e. the range of wl that contributes to chromospheric brightness)
    f_Itot_ca8 = np.sum(f_slab_ca8[7:12,:])/f_slab_ca8[7:12,:].size

    #Background extract
    #===========
    #inversion
    b_obs = b_i.dat[0,0,binit:binit+bnpix,:,ss]
    b_obs[:,0:42] = b_i.dat[0,0,binit:binit+bnpix,0:42,ss]#*amp
    b_syn = b_o.dat[0,0,binit:binit+bnpix,:,ss]
    b_syn[:,0:42] = b_o.dat[0,0,binit:binit+bnpix,0:42,ss]#*amp
    b_wav = b_i.wav
    b_n = bnpix
    b_dep = b_m.ltau[0,0,0,:]
    #slab
    bg_file = file_search(fibdir,'crispex*3950*.csav')[2*fn+1]
    bg = restore(fibdir+bg_file)
    b_slab = bg.loop_slab[:,fr,:]*1000 #intensity values in the desired frame
    b_x = bg.x_coords
    b_x_pts = bg.x_loop_pts
    b_y_pts = bg.y_loop_pts
    b_y = bg.y_coords
    b_l = bg.loop_size
    #calculating the effective intensity
    #(i.e. the range of wl that contributes to chromospheric brightness)
    b_Itot = np.sum(b_slab[6:15,:])/b_slab[6:15,:].size

    # bg in Ha
    bg_file_h = (file_search(fibdir,'crispex*6563*.csav'))[2*fn+1]
    bg_h = restore(fibdir+bg_file_h)
    b_slab_h = bg_h.loop_slab[:,fr,:]/calib_ha[0]  #intensity values in the desired frame
    #calculating the effective intensity
    #(i.e. the range of wl that contributes to chromospheric brightness)
    b_Itot_h = np.sum(b_slab_h[4:10,:])/b_slab_h[4:10,:].size

    # bg in Ca 8542
    bg_file_ca8 = (file_search(fibdir,'crispex*8542*.csav'))[2*fn+1]
    bg_ca8 = restore(fibdir+bg_file_ca8)
    b_slab_ca8 = bg_ca8.loop_slab[0,:,fr,:]/calib_ca8[0]  #intensity values in the desired frame
    #calculating the effective intensity
    #(i.e. the range of wl that contributes to chromospheric brightness)
    b_Itot_ca8 = np.sum(b_slab_ca8[7:12,:])/b_slab_ca8[7:12,:].size

        # FIBRIL and BG IDs:
    index = ((fib_file.split('_')[5]+'_'+fib_file.split('_')[6]).split('.')[0]).split('Aug')[1]
    index_bg = ((bg_file.split('_')[5]+'_'+bg_file.split('_')[6]).split('.')[0]).split('Aug')[1]
    for id in range(len(indxs)):
        if indxs[id] == index:
            print 'f, b = ', [index,index_bg], 'id = ', fn

    # Intensity difference between F and B + CALIBRATION:
    I_dif[fn] = 100.*(f_Itot - b_Itot)/b_Itot
    I_dif_h[fn] = 100.*(f_Itot_h - b_Itot_h)/b_Itot_h
    I_dif_ca8[fn] = 100.*(f_Itot_ca8 - b_Itot_ca8)/b_Itot_ca8

    # temperature of each fibril and neighboring background
    f_temp_cak = f_m.temp[0, 0, finit:finit+bnpix, approx_cak_dep_n]
    b_temp_cak = b_m.temp[0, 0, binit:binit+bnpix, approx_cak_dep_n]
    t_dif_mean[fn] = np.mean(f_temp_cak) - np.mean(b_temp_cak)
    t_dif_each[binit:binit+bnpix] = f_temp_cak - b_temp_cak    

    #break
    # Alert on the false fibril detections
    if I_dif[fn] <= 0:
        print 'False detection index = '+index
    else:
        #print 'background: ',bg_file,'    ', bnpix, '=?', len(b_y_pts)
        #containers:
        f_temp_m = np.zeros(ndep)
        f_temp_med = np.zeros(ndep)
        f_los_m = np.zeros(ndep)
        f_los_med = np.zeros(ndep)
        f_turb_m = np.zeros(ndep)
        f_turb_med = np.zeros(ndep)
        b_temp_m = np.zeros(ndep)
        b_temp_med = np.zeros(ndep)
        b_los_m = np.zeros(ndep)
        b_turb_m = np.zeros(ndep)
        b_los_med = np.zeros(ndep)
        b_turb_med = np.zeros(ndep)
        
        #chi2
        pxnl = np.min([f_n,b_n])
        chi2_f = np.zeros(pxnl)
        chi2_b = np.zeros(pxnl)
        for px in range(0,pxnl):
            chi2_f[px] = np.sum(((f_o.dat[0,0,px+finit,:,ss] - f_i.dat[0,0,px+finit,:,ss])/f_o.weights[:,ss])**2)/len(f_wav)
            chi2_b[px] = np.sum(((b_o.dat[0,0,px+binit,:,ss] - b_i.dat[0,0,px+binit,:,ss])/b_o.weights[:,ss])**2)/len(f_wav)
            chi2_fm = np.mean(chi2_f)

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

        # amplifying the cak spectral profile
        amp = 3.5 #cak amplifier
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

            f = plt.figure(figsize=(8,8))
            ax1 = plt.subplot2grid((10,6), (0,0), colspan=6, rowspan=3)
            ax2 = plt.subplot2grid((10,6), (3,0), colspan=2, rowspan=3)
            ax3 = plt.subplot2grid((10,6), (3,2), colspan=2, rowspan=3)
            ax4 = plt.subplot2grid((10,6), (3,4), colspan=2, rowspan=3)
            ax5 = plt.subplot2grid((10,6),(6,0), colspan=3,rowspan=6)
            ax6 = plt.subplot2grid((10,6),(6,3), colspan=3,rowspan=6)
            
            ax1.plot(np.mean(b_obs,axis = 0),'s',color = 'blue',markersize = 3)
            ax1.plot(np.mean(f_obs, axis = 0),'.',color = 'red',markersize = 4.5)
            ax1.plot(np.mean(b_syn,axis = 0), color='b',linewidth = 1.25)
            ax1.plot(np.mean(f_syn, axis = 0), color='r',linewidth = 1.25)
            ax1.text(2,0.75,r'$\mathrm{\overline{\chi^2_{f}}}$ = '+"%.2f" % round(np.mean(chi2_f),2) + '\n' + r'$\mathrm{\overline{\chi^2_{b}}}$ = '+"%.2f" % round(np.mean(chi2_b),2),color = 'black',horizontalalignment='left', verticalalignment='bottom', bbox={'facecolor': 'white', 'pad': 0.2,'alpha':0.25, 'boxstyle':'round'})


            for bb in range(0,b_n):
                ax2.plot(b_dep, b_m.temp[0,0,bb+binit,:].squeeze()*1e-3, 'k-',color = 'blue',alpha = alp)
                ax3.plot(b_dep, b_m.vturb[0,0,bb+binit,:].squeeze()*1.e-5, 'k-',color = 'blue',alpha = alp)
                ax4.plot(b_dep, b_m.vlos[0,0,bb+binit,:].squeeze()*1.e-5, 'k-',color = 'blue',alpha = alp)
            for dd in range(ndep):
                b_temp_m[dd] = np.mean(b_m.temp[0,0,binit:binit+bnpix,dd])
                b_temp_med[dd] = np.median(b_m.temp[0,0,binit:binit+bnpix,dd])
                b_los_m[dd] = np.mean(b_m.vlos[0,0,binit:binit+bnpix,dd])
                b_los_med[dd] = np.median(b_m.vlos[0,0,binit:binit+bnpix,dd])
                b_turb_m[dd] = np.mean(b_m.vturb[0,0,binit:binit+bnpix,dd])
                b_turb_med[dd] = np.median(b_m.vturb[0,0,binit:binit+bnpix,dd])
                
            for ff in range(0,f_n):
                f_temp_m += f_m.temp[0,0,finit+ff,:]
                ax2.plot(f_dep, f_m.temp[0,0,finit+ff,:].squeeze()*1e-3, 'k-',color = 'red',alpha = alp)
                ax3.plot(f_dep, f_m.vturb[0,0,finit+ff,:].squeeze()*1.e-5, 'k-',color = 'red',alpha = alp)
                ax4.plot(f_dep, f_m.vlos[0,0,finit+ff,:].squeeze()*1.e-5, 'k-',color = 'red',alpha = alp)
            for dd in range(ndep):
                f_temp_m[dd] = np.mean(f_m.temp[0,0,finit:finit+fnpix,dd])
                f_temp_med[dd] = np.median(f_m.temp[0,0,finit:finit+fnpix,dd])
                f_los_m[dd] = np.mean(f_m.vlos[0,0,finit:finit+fnpix,dd])
                f_los_med[dd] = np.median(f_m.vlos[0,0,finit:finit+fnpix,dd])
                f_turb_m[dd] = np.mean(f_m.vturb[0,0,finit:finit+fnpix,dd])
                f_turb_med[dd] = np.median(f_m.vturb[0,0,finit:finit+fnpix,dd])
                
            lw = 1.1 #plot linewidth
            ax2.plot(b_dep, (b_temp_m).squeeze()*1e-3, 'k-',color = 'b',linewidth = lw)
            #ax2.plot(b_dep, (b_temp_med).squeeze()*1e-3, 'k-',color = 'b',linestyle = '--',linewidth = lw)
            ax2.plot(f_dep, (f_temp_m).squeeze()*1e-3, 'k-',color = 'r',linewidth = lw)
            #ax2.plot(f_dep, (f_temp_med).squeeze()*1e-3, 'k-',color = 'r',linestyle = '--',linewidth = lw)
            for n in range(len(nodes_temp)):
                ax2.axvline(nodes_temp[n],linestyle = '--',color = 'k',linewidth = 0.2)
                ax3.plot(b_dep, (b_turb_m).squeeze()*1e-5, 'k-',color = 'b',linewidth = lw)
                #ax3.plot(b_dep, (b_turb_med).squeeze()*1e-5, 'k-',color = 'b',linestyle = '--',linewidth = lw)
                ax3.plot(f_dep, (f_turb_m).squeeze()*1e-5, 'k-',color = 'r',linewidth = lw)
                #ax3.plot(f_dep, (f_turb_med).squeeze()*1e-5, 'k-',color = 'r',linestyle = '--',linewidth = lw)
            for n in range(len(nodes_vturb)):
                ax3.axvline(nodes_vturb[n],linestyle = '--',color = 'k',linewidth = 0.2)
                ax4.plot(b_dep, (b_los_m).squeeze()*1e-5, 'k-',color = 'b',linewidth = lw)
                #ax4.plot(b_dep, (b_los_med).squeeze()*1e-5, 'k-',color = 'b',linestyle = '--',linewidth = lw)
                ax4.plot(f_dep, (f_los_m).squeeze()*1e-5, 'k-',color = 'r',linewidth = lw)
                #ax4.plot(f_dep, (f_los_med).squeeze()*1e-5, 'k-',color = 'r',linestyle = '--',linewidth = lw)
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
            ax2.set_xlabel(r'log($\mathrm{\tau}$)')
            ax2.set_ylabel(r'T [kK]')
            
            ax3.set_ylim(-0.1,10)
            ax3.set_xlim(-7.5,1)
            ax3.set_xlabel(r'log($\mathrm{\tau}$)')
            ax3.set_ylabel(r'v$_{turb}$ [km/s]')
            
            ax4.set_ylim(-7.5,7.5)
            ax4.set_xlim(-7.5,1)
            ax4.set_xlabel(r'log($\mathrm{\tau}$)')
            ax4.set_ylabel(r'v$\mathrm{_{LOS}}$ [km/s]')
            
            #cropping maps
            edge = 40
            xmin = int(np.min([f_x_pts[0],f_x_pts[-1]])-edge)
            xmax = int(np.max([f_x_pts[0],f_x_pts[-1]])+edge)
            ymin = int(np.min([f_y_pts[0],f_y_pts[-1]])-edge)
            ymax = int(np.max([f_y_pts[0],f_y_pts[-1]])+edge)
            xlabel = '[pixels]'
            ylabel = '[pixels]'
            cmap = 'gray'
            title = r'Continuum 4000 $\mathrm{\AA}$'
            titlek = r'Ca II K $\mathrm{\lambda}$-integrated'
            #cont map
            c_map_cropped = c_map[ymin:ymax,xmin:xmax]
            ax5.imshow(c_map_cropped, cmap = cmap)
            ax5.set_xlim(0, xmax-xmin)
            ax5.set_ylim(0, ymax-ymin)
            ax5.set_xlabel(xlabel)
            ax5.set_ylabel(ylabel)
            #ax5.set_title(title)
            ax5.text(0.5,0.5,r'Continuum 4000 $\mathrm{\AA}$',color = 'white',horizontalalignment='left', verticalalignment='bottom') #, bbox={'facecolor': 'black', 'pad': 10,'alpha':0.5}

            #overplotting paths on the maps
            ax5.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue')
            #ax5.plot(inv.b_x_pts[inv.pxn],inv.b_y_pts[inv.pxn], '+', color = 'blue',markersize = 15)
            #overplotting the fibril
            ax5.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'red')
            #ax5.plot(inv.f_x_pts[inv.pxn],inv.f_y_pts[inv.pxn], '+', color = 'red',markersize = 15)
            #ax5.text(0.5,120,reg, color = 'black', fontsize = 9)
            
            #cak map
            cak_map_cropped = cak_int_un[ymin:ymax,xmin:xmax]
            ax6.imshow(cak_map_cropped, cmap = cmap)
            ax6.set_xlim(0, xmax-xmin)
            ax6.set_ylim(0, ymax-ymin)
            ax6.set_xlabel(xlabel)
            ax6.set_ylabel('')
            #ax6.set_title(titlek)
            ax6.text(0.5,0.5,r'Ca II K $\mathrm{\lambda}$-integrated',color = 'white',horizontalalignment='left', verticalalignment='bottom') #, bbox={'facecolor': 'black', 'pad': 10,'alpha':0.5}
            #ax6.text(7,95,'227',color = 'white',horizontalalignment='left', verticalalignment='bottom') #, bbox={'facecolor': 'black', 'pad': 10,'alpha':0.5}
            #overplotting the fibril on the maps
            ax6.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue',alpha = 0.3,linewidth = 4)
            ax6.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'red',alpha = 0.3,linewidth = 4)
            f.set_tight_layout(True)
            plt.ion()
            plt.show()

            #f.savefig(outdir+'inv_pro_'+index+'_'+mode+'_cycle'+cyc+'.pdf')
            print 'file saved to:'+outdir+'inv_pro_'+index+'_'+mode+'_cycle'+cyc+'.pdf'
            plt.close('all')

            finit = finit + fnpix
            binit = binit + bnpix



xy = np.linspace(-100,100,100)
            

##############
# the HALPHA I plot
##############
# Excluding the False detections in Intensity difference
I_dif_true = I_dif[np.where(I_dif > 0)[0]]
I_dif_h_true = I_dif_h[np.where(I_dif > 0)[0]]

x = I_dif_true#*1e7
y = I_dif_h_true#*1e7

true_ck = I_dif_true[np.where(I_dif_h_true<=0)[0]]#*1e7
true_h = I_dif_h_true[np.where(I_dif_h_true<=0)[0]]#*1e7

plt.close('all')
fig = plt.figure(figsize = (5,5))

xmin, xmax = -2, 40
ymin, ymax = -10, 40
bins = 50
gs = GridSpec(5,5)

ax_joint = fig.add_subplot(gs[1:5,0:4])
ax_marg_x = fig.add_subplot(gs[0,0:4])
ax_marg_y = fig.add_subplot(gs[1:5,4])

ax_joint.scatter(x, y, alpha = 0.5, color = 'orangered')
ax_joint.scatter(true_ck, true_h, color = 'white', edgecolor = 'white')
ax_joint.scatter(true_ck, true_h, color = 'darkgray', edgecolor = 'dimgray',alpha = 0.7)
ax_marg_x.hist(x, color = 'orangered', rwidth = 0.95, alpha = 0.5, bins = bins)
ax_marg_y.hist(y,orientation="horizontal", color = 'orangered', alpha = 0.5, rwidth = 0.95, bins = bins)

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
ax_joint.set_xlabel(r' $\mathrm{\Delta I_{Ca K}}$ [%]', fontdict = font)
ax_joint.set_ylabel(r' $\mathrm{\Delta I_{H_{\alpha}}}$ [%]', fontdict = font)
ax_joint.tick_params(labelsize = fontsize)
ax_marg_x.tick_params(labelsize = fontsize)
ax_marg_y.tick_params(labelsize = fontsize)

# Set ax limits on marginals
ax_joint.set_xlim([xmin,xmax])
ax_joint.set_ylim([ymin,ymax])
ax_marg_y.set_ylim([ymin,ymax])
ax_marg_x.set_xlim([xmin,xmax])


pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
ax_joint.annotate('pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')

angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
ax_joint.plot(xy,xy,alpha = 0.5, color = 'black', linestyle = '--')
ax_joint.text(xmax - 10,ymax - 5,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

plt.tight_layout()
plt.subplots_adjust(left = 0.135,
                    bottom = 0.11,
                    right = 0.99,
                    top = 0.98,
                    wspace = 0.11,
                    hspace = 0.11
)
plt.show()
fig.savefig(outdir + 'I_dif_scatter_density_new.pdf', dpi = 1000)







##############
# the Ca 8542 plot
##############
# Excluding the False detections in Intensity difference
I_dif_true = I_dif[np.where(I_dif > 0)[0]]
I_dif_ca8_true = I_dif_ca8[np.where(I_dif > 0)[0]]
t_dif_true = t_dif_mean[np.where(I_dif > 0)[0]]

x = I_dif_true#*1e7
y = I_dif_ca8_true#*1e7

true_ck = I_dif_true[np.where(I_dif_ca8_true<=0)[0]]#*1e7
true_ca8 = I_dif_ca8_true[np.where(I_dif_ca8_true<=0)[0]]#*1e7

plt.close('all')
fig = plt.figure(figsize = (5,4))

xmin, xmax = -2, 40
ymin, ymax = -10, 40
bins = 50
gs = GridSpec(4,5)

ax_joint = fig.add_subplot(gs[0:4,0:4])
#ax_marg_x = fig.add_subplot(gs[0,0:4])
ax_marg_y = fig.add_subplot(gs[0:4,4])

ax_joint.scatter(x, y, alpha = 0.5, color = 'orangered')
ax_joint.scatter(true_ck, true_ca8, color = 'white', edgecolor = 'white')
ax_joint.scatter(true_ck, true_ca8, color = 'darkgray', edgecolor = 'dimgray',alpha = 0.7)
#ax_marg_x.hist(x, color = 'orangered', rwidth = 0.95, alpha = 0.5, bins = bins)
ax_marg_y.hist(y,orientation="horizontal", color = 'orangered', alpha = 0.5, rwidth = 0.95, bins = bins)

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
ax_joint.set_xlabel(r' $\mathrm{\Delta I_{Ca K}}$ [%]', fontdict = font)
ax_joint.set_ylabel(r' $\mathrm{\Delta I_{Ca IR}}$ [%]', fontdict = font)
ax_joint.tick_params(labelsize = fontsize)
ax_marg_x.tick_params(labelsize = fontsize)
ax_marg_y.tick_params(labelsize = fontsize)

# Set ax limits on marginals
ax_joint.set_xlim([xmin,xmax])
ax_joint.set_ylim([ymin,ymax])
ax_marg_y.set_ylim([ymin,ymax])
ax_marg_x.set_xlim([xmin,xmax])


pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
pearsonp = str(np.int(stat.pearsonr(x, y)[1]))
ax_joint.annotate('pearsonr = ' + pearsonr,xy=(0.05,0.92), fontsize = fontsize-1, xycoords='axes fraction')
angle = 90.-(np.arctan((ymax-ymin)/float(xmax-xmin))*360./(2*np.pi))
ax_joint.plot(xy,xy,alpha = 0.5, color = 'black', linestyle = '--')
ax_joint.text(xmax - 10,ymax - 5,'y = x', alpha = 0.5, color = 'black', rotation = angle, fontsize = 8)

plt.tight_layout()
plt.subplots_adjust(left = 0.135,
                    bottom = 0.11,
                    right = 0.99,
                    top = 0.98,
                    wspace = 0.11,
                    hspace = 0.11
)
plt.show()
fig.savefig(outdir + 'I_dif_scatter_density_ca8542.pdf', dpi = 1000)
