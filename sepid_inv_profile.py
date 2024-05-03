import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
#import spectral as s

# DECLERATIONS
indx_case = 1
#index = '06_165307'
#bg_indxs = '06_165404'
index = '06_165307'
index_bg = '06_165404'
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
alp = 0.01 #transparency of the plotted lines

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

fn = 95

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
    
# slab
fibdir = datadir+'fr'+str(fr)+'/'
fib_file = file_search(fibdir,'crispex*3950*'+index+'.csav')[0]
fib = restore(fibdir+fib_file)
f_slab = fib.loop_slab[:,fr,:]*1000  #intensity values in the desired frame
f_x = fib.x_coords
f_x_pts = fib.x_loop_pts
f_y = fib.y_coords
f_l = fib.loop_size
f_y_pts = fib.y_loop_pts
    
#fibril in Ha
fib_file_h = file_search(fibdir,'crispex*6563*'+index+'.csav')[0]
fib_h = restore(fibdir+fib_file_h)
f_slab_h = fib_h.loop_slab[:,fr,:]/calib_ha[0]  #intensity values in the desired frame
    
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
    
# slab
bg_file = file_search(fibdir,'crispex*3950*'+index_bg+'.csav')[0]
bg = restore(fibdir+bg_file)
b_slab = bg.loop_slab[:,fr,:]*1000 #intensity values in the desired frame
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

# fibril in Ha
bg_file_h = file_search(fibdir,'crispex*6563*'+index_bg+'.csav')[0]
bg_h = restore(fibdir+bg_file_h)
b_slab_h = bg_h.loop_slab[:,fr,:]/calib_ha[0]  #intensity values in the desired frame

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
                            
if(ss == 0):
    ymin = 0.05
    ymax = np.max(np.mean(f_obs, axis = 0))
elif(ss == 1 or ss == 2):
    ymin = -0.055
    ymax = 0.005
else:
    ymin = -0.03
    ymax = 0.03

#amplifying the cak spectral profile
amp = 3.5 #cak amplifier
f_obs_amp = f_i.dat[0,0,finit:finit+fnpix,:,ss]
f_obs_amp[:,0:42]  = f_i.dat[0,0,finit:finit+fnpix,0:42,ss]*amp
f_syn_amp = f_o.dat[0,0,finit:finit+fnpix,:,ss]
f_syn_amp[:,0:42]  = f_o.dat[0,0,finit:finit+fnpix,0:42,ss]*amp
b_obs_amp = b_i.dat[0,0,binit:binit+bnpix,:,ss]
b_obs_amp[:,0:42] = b_i.dat[0,0,binit:binit+bnpix,0:42,ss]*amp
b_syn_amp = b_o.dat[0,0,binit:binit+bnpix,:,ss]
b_syn_amp[:,0:42] = b_o.dat[0,0,binit:binit+bnpix,0:42,ss]*amp
            
plt.close("all")

f = plt.figure(figsize=(8,8))
ax1 = plt.subplot2grid((10,6), (0,0), colspan=6, rowspan=3)
#ax11 = plt.subplot2grid((10,6), (0,0), colspan=2, rowspan=3)
#ax12 = plt.subplot2grid((10,6), (0,0), colspan=2, rowspan=3)
#ax13 = plt.subplot2grid((10,6), (0,0), colspan=2, rowspan=3)
ax2 = plt.subplot2grid((10,6), (3,0), colspan=2, rowspan=3)
ax3 = plt.subplot2grid((10,6), (3,2), colspan=2, rowspan=3)
ax4 = plt.subplot2grid((10,6), (3,4), colspan=2, rowspan=3)
ax5 = plt.subplot2grid((10,6),(6,0), colspan=3,rowspan=6)
ax6 = plt.subplot2grid((10,6),(6,3), colspan=3,rowspan=6)
                
ax1.plot(np.mean(b_obs,axis = 0),'s',color = 'blue',markersize = 3)
ax1.plot(np.mean(f_obs, axis = 0),'.',color = 'red',markersize = 4.5)
ax1.plot(np.mean(b_syn,axis = 0), color='b',linewidth = 1.25)
ax1.plot(np.mean(f_syn, axis = 0), color='r',linewidth = 1.25)
                
                
for bb in range(0,b_n):
    ax2.plot(b_dep, b_m.temp[0,0,bb+binit,:].squeeze()*1e-3, 'k-',color = 'blue',alpha = alp)
    ax3.plot(b_dep, b_m.vturb[0,0,bb+binit,:].squeeze()*1.e-5, 'k-',color = 'blue',alpha = alp)
    ax4.plot(b_dep, b_m.vlos[0,0,bb+binit,:].squeeze()*1.e-5, 'k-',color = 'blue',alpha = alp)
for dd in range(45):
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
for dd in range(45):
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

for ii in range(gg):
    ax5.plot(b_x_pts[tp[ii]]-xmin,b_y_pts[tp[ii]]-ymin, marker = markers[ii], color = 'blue',markersize = 5)
    ax5.plot(f_x_pts[tp[ii]]-xmin,f_y_pts[tp[ii]]-ymin, marker = markers[ii], color = 'red',markersize = 5)
    ax6.plot(b_x_pts[tp[ii]]-xmin,b_y_pts[tp[ii]]-ymin, marker = markers[ii], color = 'blue',markersize = 5,markerfacecolor = 'none')
    ax6.plot(f_x_pts[tp[ii]]-xmin,f_y_pts[tp[ii]]-ymin, marker = markers[ii], color = 'red',markersize = 5,markerfacecolor = 'none')

f.set_tight_layout(True)
plt.ion()
plt.show()

outname = outdir+'test_inv_pro_'+index+'_'+mode+'_cycle'+cyc+'.pdf'
f.savefig(outname)
print 'file saved to:'+outname
