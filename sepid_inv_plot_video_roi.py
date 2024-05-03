import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec

#import spectral as s

#DECLERATIONS
index = '06_165307'
index_bg = '06_165404'

mode = 'roinew_fullstokes_smoothed'
cyc = '2'
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

font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }


# settings board!
reg = 'n$_{T}$ = '+" ".join(str('%.1f'%(x)) for x in nodes_temp)+' \nn$_{v_{LOS}}$ = '+" ".join(str('%.1f'%(x)) for x in nodes_vlos)+' \nv$_{turb}$ = '+" ".join(str('%.1f'%(x)) for x in nodes_vturb)+'  \nreg$_{type}$ = '+" ".join(str(int(x)) for x in regularization_type)+'        reg = '+str('%.1f'%(regularize))+'       reg$_{weight}$ = '+" ".join(str('%.1f'%(x)) for x in regularization_weights)

alp = 0.01 #transparency of the plotted lines
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
outdir = '/home/seki2695/OUTPUT/inv/'
invdir = '/home/seki2695/INV/stic/roi/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
fn = 28 #frame no.
ss = 0 # stokes param index

# data extract
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
fn = 28 #frame no.
file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
file_ca8 = file_search(datadir,'crispex*'+pref[1]+'*.fcube')
file_cak = file_search(datadir,'crispex*'+pref[2]+'*.fcube')
file_ha = file_search(datadir,'crispex*'+pref[3]+'*.fcube')
cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])[fn]
cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])[fn]
cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])[fn]
cube_ha = lp_read(datadir+file_ha[0],datadir+file_ha[1])[fn]
c_map = cube_cak[0,-1,:,:]
cak_int_un = (lp_read_scan(datadir+'cak_int_un.fcube'))[fn,:,:]

# inversion extract
f_i = sp.profile(invdir+'roinew_observed_'+index+'.nc') # f_i.dat = (nt, ny, nx, nw, ns)
f_o = sp.profile(invdir+mode+'_synthetic_cycle'+cyc+'_'+index+'.nc')
f_m = sp.model(invdir+mode+'_atmosout_cycle'+cyc+'_'+index+'.nc')

temp = f_m.temp[0,:,:,:].squeeze()*1e-3
dep = f_m.ltau[0,0,0,:]
cont = f_i.dat[0,:,:,42,0]
wave = f_i.wav
vlos = f_m.vlos[0,:,:,:].squeeze()*1e-5
vturb = f_m.vturb[0,:,:,:].squeeze()*1e-5
dy = len(f_m.temp[0,:,0,0])
dx = len(f_m.temp[0,0,:,0])
dw = len(wave)

# chi2 map
chi2 = np.zeros((dy,dx,dw), dtype = 'float')
for xx in range(dx):
    for yy in range(dy):
        chi2[yy,xx,:] = ((f_o.dat[0,yy,xx,:,ss] - f_i.dat[0,yy,xx,:,ss])/f_o.weights[:,ss])**2
chi2_cak = np.sum(chi2[:,:,0:52],axis=2)/len(chi2[0,0,0:52])
chi2_ca8 = np.sum(chi2[:,:,53:163], axis=2)/len(chi2[0,0,53:163])
chi2_fe = np.sum(chi2[:,:,164:], axis=2)/len(chi2[0,0,164:])

# stokes v map
fe1_b0, fe1_bn, fe1_r0, fe1_rn = 0, 4, 5, 9 # blue and red wings of 1st line
fe2_b0, fe2_bn, fe2_r0, fe2_rn = 10, 13, 14, 16 # 2nd line
fe_cont_idx = 9

CP_map = (np.sum(cube_fe[3, fe1_b0: fe1_bn, :, :], axis = 0)\
         + np.sum(cube_fe[3,fe2_b0:fe2_bn, :, :], axis = 0)\
         - np. sum(cube_fe[3, fe1_r0:fe1_rn, :, :], axis = 0)\
         - np.sum(cube_fe[3,fe2_r0:fe2_rn, :, :], axis = 0))/cube_fe[0,fe_cont_idx,:,:]

CPch_map = (np.sum(cube_ca8[3,0:10,:,:],axis=0)-np.sum(cube_ca8[3,11:20,:,:],axis = 0))/cube_ca8[0,0,:,:]

# slab extract
# Fibril
fibdir = datadir+'fr'+str(fn)+'/'
fib_file = file_search(fibdir,'crispex*3950*'+index+'*.csav')
fib = restore(fibdir+fib_file[0])
f_slab = fib.loop_slab  #intensity values in the desired frame
f_x = fib.x_coords
f_x_pts = fib.x_loop_pts
f_y = fib.y_coords
f_l = fib.loop_size
f_y_pts = fib.y_loop_pts

#BackGround
bg_file = file_search(fibdir,'crispex*3950*'+index_bg+'*.csav')
bg = restore(fibdir+bg_file[0])
b_slab = bg.loop_slab #intensity values in the desired frame
b_x = bg.x_coords
b_x_pts = bg.x_loop_pts
b_y_pts = bg.y_loop_pts
b_y = bg.y_coords
b_l = bg.loop_size

# roi coords
edge = 40
xmin = int(np.min([f_x_pts[0],f_x_pts[-1]])-edge)
xmax = int(np.max([f_x_pts[0],f_x_pts[-1]])+edge)
ymin = int(np.min([f_y_pts[0],f_y_pts[-1]])-edge)
ymax = int(np.max([f_y_pts[0],f_y_pts[-1]])+edge)

# preparing for the video
#fr_i, fr_n = 5, len(dep)-5
fr_i, fr_n = 7, len(dep)-7

for dd in range(fr_i, fr_n):
#for dd in range(1):
  #  dd = 21
    plt.close('all')

    f = plt.figure(figsize = [5.,6])
    gs = gridspec.GridSpec(3, 1) # grid scale of the 1st col.
    #plt.style.use('v2.0')
    aspect =1# float(f_n_valid)/float(len(f_depRF_clipped))/2.5 #len(f_depRF_clipped)/f_n_valid # panels aspect ratio (y/x)
    gs.update(left=0.075,
              right=0.88,
              #wspace=0.05,
              bottom=0.08,
              top=0.99,
              hspace = 0.0,
    )

    xticks = np.round(np.linspace(0,280,280/40+1))
    yticks = np.round(np.linspace(0,120,120/30+1))
    ax1 = plt.subplot(gs[0,0],adjustable = 'box')#,aspect = 'equal')
    ax1.set_xticks([])
    ax1.set_ylim(0,dy)
    ax1.set_yticks(yticks)
    ax1.set_ylabel(r'[pixels]', fontdict=font)
    ax2 = plt.subplot(gs[1,0],adjustable = 'box')#,aspect = 'equal')
    ax2.set_xticks([])
    ax2.set_ylim(0,dy)
    ax2.set_yticks(yticks)
    ax2.set_ylabel(r'[pixels]', fontdict=font)
    ax3 = plt.subplot(gs[2,0],adjustable = 'box')#,aspect = 'equal')
    ax3.set_xlabel(r'[pixels]', fontdict=font)
    ax3.set_ylim(0,dy)
    ax3.set_yticks(yticks)
    ax3.set_xticks(xticks)
    ax3.set_ylabel(r'[pixels]', fontdict=font)
    
    # temp colorbar axis
    axin1 = inset_axes(ax1,
                       width="3%",  # width = 10% of parent_bbox width
                       height="90%",  # height : 50%
                       loc='center left',
                       bbox_to_anchor=(1.021, 0., 1, 1),
                       bbox_transform=ax1.transAxes,
                       borderpad=0,
    )
    
    # vlos colorbar axis
    axin2 = inset_axes(ax2, #axis
                       width="3%",  # width = 10% of parent_bbox width
                       height="90%",  # height : 50%
                       loc='center left',
                       bbox_to_anchor=(1.02, 0., 1, 1),
                       bbox_transform=ax2.transAxes,
                       borderpad=0,
    )
    
    # vmt colorbar axis
    axin3 = inset_axes(ax3, #axis
                       width="3%",  # width = 10% of parent_bbox width
                       height="90%",  # height : 50%
                       loc='center left',
                       bbox_to_anchor=(1.02, 0., 1, 1),
                       bbox_transform=ax3.transAxes,
                       borderpad=0,
    )
    
    # temperature panel
    panel_temp = ax1.imshow(temp[:,:,dd],cmap = 'inferno', vmax = 7)
    # overplotting the fibril
    #ax1.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'black',alpha = 0.75,linewidth = 0.5, linestyle = '--')
    ax1.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'white',alpha = 0.75,linewidth = 0.75, linestyle = '--')

    #temp_min, temp_max = np.min(temp[:,:,dd]), np.max(temp[:,:,dd])
    #ct_ticks = np.linspace(np.round(temp_min),np.round(temp_max),np.abs(np.round(temp_min)-np.round(temp_max))+1,dtype='uint')
    ct = plt.colorbar(panel_temp, cax=axin1, orientation="vertical")
    ct.set_label(r'T (kK)', fontdict = font)
    #ct.ax.set_yticklabels(map(str,ct_ticks[:-1])+['>'+str(ct_ticks[-1])], fontdict=font)
    
    # vlos panel
    norm_vlos = np.max(np.abs(vlos[:,:,dd]))
    panel_vlos = ax2.imshow(vlos[:,:,dd], cmap = 'bwr', vmin = -norm_vlos, vmax = norm_vlos)
    #ax2.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'black',alpha = 0.75,linewidth = 0.5, linestyle = '--')
    ax2.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'black',alpha = 0.75,linewidth = 0.75, linestyle = '--')
    cv = plt.colorbar(panel_vlos, cax=axin2, orientation="vertical")
    cv.set_label(r'v$_{\rm LOS}$ [km/s]', fontdict=font)
    
    # vturb panel
    panel_vturb = ax3.imshow(vturb[:,:,dd], cmap = 'bone')
    #ax3.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'black',alpha = 0.75,linewidth = 0.5, linestyle = '--')
    ax3.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'white',alpha = 0.75,linewidth = 0.75, linestyle = '--')
    cvmt = plt.colorbar(panel_vturb, cax=axin3, orientation="vertical")
    cvmt.set_label(r'v$_{\rm turb}$ [km/s]', fontdict = font)
    ax1.text(5,dy-17,r'log($\tau_{500}$) = '+str(np.round(dep[dd],decimals = 2)),color = 'black',horizontalalignment='left', verticalalignment='bottom', bbox={'facecolor': 'white', 'pad': 3,'alpha':0.75})
    plt.show()

    break
    # to make the video staring from the photosphere
    #f.savefig(outdir+'video/'+'inv_'+index+'_'+mode+'_cycle'+cyc+'fr'+str(fr_n-(dd)-1).zfill(4)+'.png',dpi = 1000)

# creating the video
#os.system("ffmpeg -r 2 -i "+outdir+'video/'+'inv_'+index+'_'+mode+'_cycle'+cyc+'fr'+"%04d.png -vcodec mpeg4 -y "+outdir+'video/inv_'+index+'_'+mode+'_fullstokes_cycle'+cyc+"_new.mp4")

'''
# maps
plt.close('all')

f = plt.figure(figsize = [7.5,4.05])
gs = gridspec.GridSpec(2, 2) # grid scale of the 1st col.
aspect =1# float(f_n_valid)/float(len(f_depRF_clipped))/2.5 #len(f_depRF_clipped)/f_n_valid # panels aspect ratio (y/x)
gs.update(left=0.08,
          right=0.99,
          wspace=0.0,
          bottom=0.1,
          top=0.99,
          hspace = 0.0,
)

xticks = np.round(np.linspace(0,270,270/30+1))
yticks = np.round(np.linspace(0,120,120/30+1))

ax1 = plt.subplot(gs[0,0],adjustable = 'box')#,aspect = 'equal')
ax1.set_xticks([])
ax1.set_ylim(0,dy)
ax1.set_yticks(yticks)
ax1.set_ylabel(r'[pixels]', fontdict=font)
ax1.set_xlim(0, xmax-xmin)
ax1.set_ylim(0, ymax-ymin)

ax2 = plt.subplot(gs[0,1],adjustable = 'box',aspect = 'equal')
ax2.set_xticks([])
ax2.set_ylim(0,dy)
ax2.set_yticks([])
#ax2.set_ylabel(r'[pixels]', fontdict=font)
ax2.set_xlim(0, xmax-xmin)
ax2.set_ylim(0, ymax-ymin)

ax3 = plt.subplot(gs[1,0],adjustable = 'box',aspect = 'equal')
ax3.set_xlabel(r'[pixels]', fontdict=font)
ax3.set_ylim(0,dy)
ax3.set_yticks(yticks)
ax3.set_xticks(xticks)
ax3.set_ylabel(r'[pixels]', fontdict=font)
ax3.set_xlim(0, xmax-xmin)
ax3.set_ylim(0, ymax-ymin)

ax4 = plt.subplot(gs[1,1],adjustable = 'box',aspect = 'equal')
ax4.set_xlabel(r'[pixels]', fontdict=font)
ax4.set_ylim(0,dy)
ax4.set_yticks([])
ax4.set_xticks(xticks)
#ax4.set_ylabel(r'[pixels]', fontdict=font)
ax4.set_xlim(0, xmax-xmin)
ax4.set_ylim(0, ymax-ymin)


cmap = 'gray'
    
# ax1 colorbar axis
#axin1 = inset_axes(ax1,
#                   width="3%",  # width = 10% of parent_bbox width
#                   height="90%",  # height : 50%
#                   loc='center left',
#                   bbox_to_anchor=(1.021, 0., 1, 1),
#                   bbox_transform=ax1.transAxes,
#                   borderpad=0,
#)
    
# ax2 colorbar axis
#axin2 = inset_axes(ax2, #axis
#                   width="3%",  # width = 10% of parent_bbox width
#                   height="90%",  # height : 50%
#                   loc='center left',
#                   bbox_to_anchor=(1.02, 0., 1, 1),
#                   bbox_transform=ax2.transAxes,
#                   borderpad=0,
#)
    
# ax3 colorbar axis
#axin3 = inset_axes(ax3, #axis
#                   width="3%",  # width = 10% of parent_bbox width
#                   height="90%",  # height : 50%
#                   loc='center left',
#                   bbox_to_anchor=(1.02, 0., 1, 1),
#                   bbox_transform=ax3.transAxes,
#                   borderpad=0,
#)


# cont map
c_map_cropped = c_map[ymin:ymax,xmin:xmax]
ax1.imshow(c_map_cropped, cmap = cmap)
ax1.text(0.5,0.5,r'Continuum 4000 $\mathrm{\AA}$',color = 'white',horizontalalignment='left', verticalalignment='bottom') #, bbox={'facecolor': 'black', 'pad': 10,'alpha':0.5}
# overplotting F and B
ax1.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue')
ax1.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'red')

# cak map
cak_map_cropped = cak_int_un[ymin:ymax,xmin:xmax]
ax2.imshow(cak_map_cropped, cmap = cmap)
ax2.text(0.5,0.5,r'Ca II K $\mathrm{\lambda}$-integrated',color = 'white',horizontalalignment='left', verticalalignment='bottom') 
# overplotting F and B
ax2.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue',alpha = 0.15,linewidth = 2)
ax2.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'red',alpha = 0.15,linewidth = 2)

# CP map at photosphere
cp_map_cropped = CP_map[ymin:ymax,xmin:xmax]
ax3.imshow(cp_map_cropped, cmap = cmap)
ax3.text(0.5,0.5,r'CP at photosphere',color = 'black',horizontalalignment='left', verticalalignment='bottom') 
# overplotting F and B
ax3.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue',alpha = 0.15,linewidth = 2)
ax3.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'red',alpha = 0.15,linewidth = 2)

# CP map at chromosphere
cpch_map_cropped = CPch_map[ymin:ymax,xmin:xmax]
ax4.imshow(cpch_map_cropped, cmap = cmap)
ax4.text(0.5,0.5,r'CP at Chromosphere',color = 'black',horizontalalignment='left', verticalalignment='bottom') 
# overplotting F and B
ax4.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue',alpha = 0.15,linewidth = 2)
ax4.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'red',alpha = 0.15,linewidth = 2)

#f.set_tight_layout(True)
plt.ion()
plt.show()

f.savefig(outdir+'inv_roi_'+index+'_'+mode+'_cycle'+cyc+'.pdf', quality = 100)
print 'file saved to:'+outdir+'inv_roi_'+index+'_'+mode+'_cycle'+cyc+'.pdf'
'''
'''
plt.close('all')
f = plt.figure(figsize = [5,6])
gs = gridspec.GridSpec(3, 1) # grid scale of the 1st col.
plt.style.use('dark_background')
font = {'family': 'sans-serif',
        'color':  'white',
        'weight': 'normal',
        'size': 8.,
}
aspect =1# float(f_n_valid)/float(len(f_depRF_clipped))/2.5 #len(f_depRF_clipped)/f_n_valid # panels aspect ratio (y/x)
gs.update(left=0.075,
          right=0.88,
          #wspace=0.05,
          bottom=0.08,
          top=0.99,
          hspace = 0.0,
)

cmap = 'gnuplot2'
xticks = np.round(np.linspace(0,280,270/40+1))
yticks = np.round(np.linspace(0,120,120/30+1))

ax1 = plt.subplot(gs[0,0],adjustable = 'box')#,aspect = 'equal')
ax1.set_xticks([])
ax1.set_ylim(0,dy)
ax1.set_yticks(yticks)
ax1.set_ylabel(r'[pixels]', fontdict=font)
ax1.set_xlim(0, xmax-xmin)
ax1.set_ylim(0, ymax-ymin)

ax2 = plt.subplot(gs[1,0],adjustable = 'box',aspect = 'equal')
ax2.set_xticks([])
ax2.set_ylim(0,dy)
ax2.set_yticks(yticks)
ax2.set_ylabel(r'[pixels]', fontdict=font)
ax2.set_xlim(0, xmax-xmin)
ax2.set_ylim(0, ymax-ymin)

ax3 = plt.subplot(gs[2,0],adjustable = 'box',aspect = 'equal')
ax3.set_xlabel(r'[pixels]', fontdict=font)
ax3.set_ylim(0,dy)
ax3.set_yticks(yticks)
ax3.set_xticks(xticks)
ax3.set_ylabel(r'[pixels]', fontdict=font)
ax3.set_xlim(0, xmax-xmin)
ax3.set_ylim(0, ymax-ymin)

vmin = np.min([chi2_cak,chi2_ca8, chi2_fe])
vmax = np.max([chi2_cak,chi2_ca8, chi2_fe])
# cak chi2
cak_chi2 = ax1.imshow(chi2_cak, cmap = cmap, vmin = vmin, vmax = vmax)
ax1.text(4.,0.5,r'Ca II K',color = 'white',horizontalalignment='left', verticalalignment='bottom') #, bbox={'facecolor': 'black', 'pad': 10,'alpha':0.5}
# overplotting F and B
#ax1.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue',alpha = 0.15,linewidth = 2)
ax1.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'white',alpha = 0.35,linewidth = 1.2, linestyle = '--')

# ca8 chi2
ca8_chi2 = ax2.imshow(chi2_ca8, cmap = cmap, vmin = vmin, vmax = vmax)
ax2.text(4,0.5,r'Ca 8542',color = 'white',horizontalalignment='left', verticalalignment='bottom') 
# overplotting F and B
#ax2.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue',alpha = 0.15,linewidth = 2)
ax2.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'white',alpha = 0.35,linewidth = 1.2, linestyle = '--')

# fe chi2
fe_chi2 = ax3.imshow(chi2_fe, cmap = cmap, vmin = vmin, vmax = vmax)
ax3.text(4,0.5,r'Fe I',color = 'white',horizontalalignment='left', verticalalignment='bottom') 
# overplotting F and B
#ax3.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'blue',alpha = 0.15,linewidth = 2)
ax3.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'white',alpha = 0.35,linewidth = 1.2, linestyle = '--')

# ax1 colorbar axis
axin1 = inset_axes(ax3,
                   width="3%",  # width = 10% of parent_bbox width
                   height="90%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.021, 0., 1, 3.335),
                   bbox_transform=ax3.transAxes,
                   borderpad=0,
)

ct = plt.colorbar(fe_chi2, cax=axin1, orientation="vertical")
ct.set_label(r'$\chi^{2}$', fontdict = font)

#f.set_tight_layout(True)
plt.ion()
plt.show()

f.savefig(outdir+'inv_roi_chi2_'+index+'_'+mode+'_cycle'+cyc+'.pdf', quality = 100)
print 'file saved to:'+outdir+'inv_roi_chi2_'+index+'_'+mode+'_cycle'+cyc+'.pdf'
'''
