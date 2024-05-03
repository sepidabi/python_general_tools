import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
#import spectral as s

#DECLERATIONS
index = '06_165307'
index_bg = '06_165404'
mod = 'XXIX'
reg = 'n$_{T}$ = -6.8, -6.0, -5.2, -4.3, -3.4, -2.5, -1.6, -0.08, 0.08 \nn$_{v_{LOS}}$ = -6.8, -4.5, -2.1, 0.2 \nv$_{turb}$ = 0.1 m/s  \nreg$_{type}$ = 1,2,3,1,1,0,1        reg = 1.0       reg$_{weight}$ = 1, 1, 1, 1, 1, 1, 1'

cyc = '3'
alp = 1. #transparency of the plotted lines
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
outdir = '/home/seki2695/OUTPUT/inv/video/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
fn = 28 #frame no.
ss = 0 # stokes param index
amp = 1.
#xx = 100  # set to zero to choose the first pixel
#yy = 0     #stokes I values in the dat array


#Extracting inversion results
f_i = sp.profile('fb_test/f'+mod+'_observed2_'+index+'.nc')
f_o = sp.profile('fb_test/f'+mod+'_synthetic_cycle'+cyc+'_'+index+'.nc')
f_m = sp.model('fb_test/f'+mod+'_atmosout_cycle'+cyc+'_'+index+'.nc')
b_i = sp.profile('bg/b'+mod+'_observed2_'+index+'.nc')
b_o = sp.profile('bg/b'+mod+'_synthetic_cycle'+cyc+'_'+index+'.nc')
b_m = sp.model('bg/b'+mod+'_atmosout_cycle'+cyc+'_'+index+'.nc')

#Fibril extract
##========
#inversion
f_obs = f_i.dat[0,0,:,:,ss]
f_syn = f_o.dat[0,0,:,:,ss]
f_wav = f_i.wav
f_n = f_i.dat.shape[2]
f_dep = f_m.ltau[0,0,0,:]
#slab
fibdir = datadir+'fr'+str(fn)+'/'
fib_file = file_search(fibdir,'crispex*3950*'+index+'*.csav')
fib = restore(fibdir+fib_file[0])
f_slab = fib.loop_slab  #intensity values in the desired frame
f_x = fib.x_coords
f_x_pts = fib.x_loop_pts
f_y = fib.y_coords
f_l = fib.loop_size
f_y_pts = fib.y_loop_pts

#Background extract
#===========
#inversion
b_obs = b_i.dat[0,0,:,:,ss]
b_syn = b_o.dat[0,0,:,:,ss]
b_wav = b_i.wav
b_n = b_i.dat.shape[2]
b_dep = b_m.ltau[0,0,0,:]
#slab
bg_file = file_search(fibdir,'crispex*3950*'+index_bg+'*.csav')
bg = restore(fibdir+bg_file[0])
b_slab = bg.loop_slab #intensity values in the desired frame
b_x = bg.x_coords
b_x_pts = bg.x_loop_pts
b_y_pts = bg.y_loop_pts
b_y = bg.y_coords
b_l = bg.loop_size

#containers:
khali = np.zeros(45)
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

#mean & median profiles
depn = len(f_dep) #number of depths
for dd in range(depn):
    b_temp_m[dd] = np.mean(b_m.temp[0,0,:,dd])
    b_temp_med[dd] = np.median(b_m.temp[0,0,:,dd])
    b_los_m[dd] = np.mean(b_m.vlos[0,0,:,dd])
    b_los_med[dd] = np.median(b_m.vlos[0,0,:,dd])
    b_turb_m[dd] = np.mean(b_m.vturb[0,0,:,dd])
    b_turb_med[dd] = np.median(b_m.vturb[0,0,:,dd])

for dd in range(depn):
    f_temp_m[dd] = np.mean(f_m.temp[0,0,:,dd])
    f_temp_med[dd] = np.median(f_m.temp[0,0,:,dd])
    f_los_m[dd] = np.mean(f_m.vlos[0,0,:,dd])
    f_los_med[dd] = np.median(f_m.vlos[0,0,:,dd])
    f_turb_m[dd] = np.mean(f_m.vturb[0,0,:,dd])
    f_turb_med[dd] = np.median(f_m.vturb[0,0,:,dd])

#maps
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

#chi2
pxnl = np.min([f_n,b_n])
chi2_f = np.zeros(pxnl)
chi2_b = np.zeros(pxnl)
for px in range(0,pxnl):
    chi2_f[px] = np.sum(((f_o.dat[0,0,px,:,ss] - f_i.dat[0,0,px,:,ss])/f_o.weights[:,ss])**2)/len(f_wav)
    chi2_b[px] = np.sum(((b_o.dat[0,0,px,:,ss] - b_i.dat[0,0,px,:,ss])/b_o.weights[:,ss])**2)/len(f_wav)
chi2_fm = np.mean(chi2_f)

#amplifying the cak spectral profile
amp = 3.5 #cak amplifier
f_obs_amp = f_i.dat[0,0,:,:,ss]
f_obs_amp[:,0:42]  = f_i.dat[0,0,:,0:42,ss]*amp
f_syn_amp = f_o.dat[0,0,:,:,ss]
f_syn_amp[:,0:42]  = f_o.dat[0,0,:,0:42,ss]*amp
b_obs_amp = b_i.dat[0,0,:,:,ss]
b_obs_amp[:,0:42] = b_i.dat[0,0,:,0:42,ss]*amp
b_syn_amp = b_o.dat[0,0,:,:,ss]
b_syn_amp[:,0:42] = b_o.dat[0,0,:,0:42,ss]*amp

#PLOTTING Profiles pixel by pixel
#===================
pxnl = np.min([f_n,b_n])
#pxnl = 5
for pxn in range(0,pxnl):
    plt.close("all")

    f = plt.figure(figsize=(8,8))
    ax1 = plt.subplot2grid((10,6), (0,0), colspan=6, rowspan=3)
    ax2 = plt.subplot2grid((10,6), (3,0), colspan=2, rowspan=3)
    ax3 = plt.subplot2grid((10,6), (3,2), colspan=2, rowspan=3)
    ax4 = plt.subplot2grid((10,6), (3,4), colspan=2, rowspan=3)
    ax5 = plt.subplot2grid((10,6),(6,0), colspan=3,rowspan=6)
    ax6 = plt.subplot2grid((10,6),(6,3), colspan=3,rowspan=6)

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
    
    #spectral profile plot
    ax1.plot(b_obs_amp[pxn,:],'s',color = 'blue',markersize = 3)
    ax1.plot(f_obs_amp[pxn,:],'.',color = 'red',markersize = 4.5)
    ax1.plot(b_syn_amp[pxn,:], color='b',linewidth = 1.25)
    ax1.plot(f_syn_amp[pxn,:], color='r',linewidth = 1.25)
    ax1.text(2,0.83,r'$\mathrm{\chi^2_{f}}$ = '+"%.2f" % round(chi2_f[pxn],2) + '\n' + r'$\mathrm{\chi^2_{b}}$ = '+"%.2f" % round(chi2_b[pxn],2),color = 'black',horizontalalignment='left', verticalalignment='bottom', bbox={'facecolor': 'white', 'pad': 0.2,'alpha':0.5, 'boxstyle':'round'})
    #break

    lw = 1. #plot linewidth
    ax2.plot(b_dep, (b_temp_m).squeeze()*1e-3, 'k-',color = 'b',linewidth = lw,alpha = 0.25)
    ax2.plot(b_dep, (b_temp_med).squeeze()*1e-3, 'k-',color = 'b',linestyle = '--',linewidth = lw,alpha = 0.25)
    ax2.plot(f_dep, (f_temp_m).squeeze()*1e-3, 'k-',color = 'r',linewidth = lw,alpha = 0.25)
    ax2.plot(f_dep, (f_temp_med).squeeze()*1e-3, 'k-',color = 'r',linestyle = '--',linewidth = lw,alpha = 0.25)
    ax3.plot(b_dep, (b_turb_m).squeeze()*1e-5, 'k-',color = 'b',linewidth = lw,alpha = 0.25)
    ax3.plot(b_dep, (b_turb_med).squeeze()*1e-5, 'k-',color = 'b',linestyle = '--',linewidth = lw,alpha = 0.25)
    ax3.plot(f_dep, (f_turb_m).squeeze()*1e-5, 'k-',color = 'r',linewidth = lw,alpha = 0.25)
    ax3.plot(f_dep, (f_turb_med).squeeze()*1e-5, 'k-',color = 'r',linestyle = '--',linewidth = lw,alpha = 0.25)
    ax4.plot(b_dep, (b_los_m).squeeze()*1e-5, 'k-',color = 'b',linewidth = lw,alpha = 0.25)
    ax4.plot(b_dep, (b_los_med).squeeze()*1e-5, 'k-',color = 'b',linestyle = '--',linewidth = lw,alpha = 0.25)
    ax4.plot(f_dep, (f_los_m).squeeze()*1e-5, 'k-',color = 'r',linewidth = lw,alpha = 0.25)
    ax4.plot(f_dep, (f_los_med).squeeze()*1e-5, 'k-',color = 'r',linestyle = '--',linewidth = lw,alpha = 0.25)

    ax2.plot(b_dep, b_m.temp[0,0,pxn,:].squeeze()*1e-3, 'k-',color = 'blue',alpha = alp)
    ax2.plot(f_dep, f_m.temp[0,0,pxn,:].squeeze()*1e-3, 'k-',color = 'red',alpha = alp)
    ax3.plot(b_dep, b_m.vturb[0,0,pxn,:].squeeze()*1.e-5, 'k-',color = 'blue',alpha = alp)
    ax3.plot(f_dep, f_m.vturb[0,0,pxn,:].squeeze()*1.e-5, 'k-',color = 'red',alpha = alp)
    ax4.plot(b_dep, b_m.vlos[0,0,pxn,:].squeeze()*1.e-5, 'k-',color = 'blue',alpha = alp)
    ax4.plot(f_dep, f_m.vlos[0,0,pxn,:].squeeze()*1.e-5, 'k-',color = 'red',alpha = alp)    

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
    xmin = int(np.min([f_x_pts[0],f_x_pts[-1]])-20)
    xmax = int(np.max([f_x_pts[0],f_x_pts[-1]])+20)
    ymin = int(np.min([f_y_pts[0],f_y_pts[-1]])-20)
    ymax = int(np.max([f_y_pts[0],f_y_pts[-1]])+20)
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
    ax5.text(7,95,pxn,color = 'red',horizontalalignment='left', verticalalignment='bottom') #, bbox={'facecolor': 'black', 'pad': 10,'alpha':0.5}
    #overplotting the background
    ax5.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue')
    ax5.plot(b_x_pts[pxn]-xmin,b_y_pts[pxn]-ymin,'s',color = 'blue',markersize = 5.,markerfacecolor = 'none')

    #overplotting the fibril
    ax5.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'red')
    ax5.plot(f_x_pts[pxn]-xmin,f_y_pts[pxn]-ymin,'s',color = 'red',markersize = 5.,markerfacecolor = 'none')
    ax5.text(0.5,120,reg + "%.2f" % round(np.mean(chi2_f),2), color = 'black', fontsize = 9)
    
    
    #cak map
    cak_map_cropped = cak_int_un[ymin:ymax,xmin:xmax]
    ax6.imshow(cak_map_cropped, cmap = cmap)
    ax6.set_xlim(0, xmax-xmin)
    ax6.set_ylim(0, ymax-ymin)
    ax6.set_xlabel(xlabel)
    ax6.set_ylabel('')
    #ax6.set_title(titlek)
    ax6.text(0.5,0.5,r'Ca II K $\mathrm{\lambda}$-integrated',color = 'white',horizontalalignment='left', verticalalignment='bottom') #, bbox={'facecolor': 'black', 'pad': 10,'alpha':0.5}
    ax6.text(7,95,pxn,color = 'red',horizontalalignment='left', verticalalignment='bottom') #, bbox={'facecolor': 'black', 'pad': 10,'alpha':0.5}
    
    #overplotting the fibril & bg
    ax6.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue',alpha = 0.3,linewidth = 4)
    ax6.plot(b_x_pts[pxn]-xmin,b_y_pts[pxn]-ymin,'s',color = 'blue',markersize = 5.,markerfacecolor = 'none')    
    ax6.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'red',alpha = 0.3,linewidth = 4)
    ax6.plot(f_x_pts[pxn]-xmin,f_y_pts[pxn]-ymin,'s',color = 'red',markersize = 5.,markerfacecolor = 'none')
    
    f.set_tight_layout(True)
    plt.ion()
    plt.show()

    f.savefig(outdir+'inv_pro_'+index+'_'+mod+'_cycle'+cyc+'_fr'+str(pxn).zfill(4)+'.png',dpi = 1000)

os.system("ffmpeg -r 4 -i "+outdir+'inv_pro_'+index+'_'+mod+'_cycle'+cyc+'_fr'+"%04d.png -vcodec mpeg4 -y "+outdir+'inv_pro_'+index+'_'+mod+'_cycle'+cyc+".mp4")
