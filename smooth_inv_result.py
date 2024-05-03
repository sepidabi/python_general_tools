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

mode = 'roinew'
cyc = '1'
# Same settings as XXVI -> for actual inversion
nodes_temp = -6.88333, -5.96667, -5.23333, -4.31667, -3.4, -2.48333, -1.56667, -0.833333, 0.083333
nodes_vlos =  -6.8, -4.46666667, -2.13333333, 0.2
nodes_vturb = 0#-6.0, -4, -1.9, -0.2
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

alp = 0.01 #transparency of the plotted lines
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
outdir = '/home/seki2695/OUTPUT/inv/'
invdir = '/scratch/sepid/stic/fov/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
fn = 28 #frame no.
ss = 0 # stokes param index

# settings board!
#reg = 'n$_{T}$ = '+" ".join(str('%.1f'%(x)) for x in nodes_temp)+' \nn$_{v_{LOS}}$ = '+" ".join(str('%.1f'%(x)) for x in nodes_vlos)+' \nv$_{turb}$ = '+" ".join(str('%.1f'%(x)) for x in nodes_vturb)+'  \nreg$_{type}$ = '+" ".join(str(int(x)) for x in regularization_type)+'        reg = '+str('%.1f'%(regularize))+'       reg$_{weight}$ = '+" ".join(str('%.1f'%(x)) for x in regularization_weights)

# data extract
#pref = ['6302','8542','3950','6563']
#pref_name = ['fe','ca8','cak','ha']
#fn = 28 #frame no.
#file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
#file_ca8 = file_search(datadir,'crispex*'+pref[1]+'*.fcube')
#file_cak = file_search(datadir,'crispex*'+pref[2]+'*.fcube')
#file_ha = file_search(datadir,'crispex*'+pref[3]+'*.fcube')
#cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])[fn]
#cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])[fn]
#cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])[fn]
#cube_ha = lp_read(datadir+file_ha[0],datadir+file_ha[1])[fn]
#c_map = cube_cak[0,-1,:,:]
#cak_int_un = (lp_read_scan(datadir+'cak_int_un.fcube'))[fn,:,:]

# slab extract
# Fibril
#fibdir = datadir+'fr'+str(fn)+'/'
#fib_file = file_search(fibdir,'crispex*3950*'+index+'*.csav')
#fib = restore(fibdir+fib_file[0])
#f_slab = fib.loop_slab  #intensity values in the desired frame
#f_x = fib.x_coords
#f_x_pts = fib.x_loop_pts
#f_y = fib.y_coords
#f_l = fib.loop_size
#f_y_pts = fib.y_loop_pts

#BackGround
#bg_file = file_search(fibdir,'crispex*3950*'+index_bg+'*.csav')
#bg = restore(fibdir+bg_file[0])
#b_slab = bg.loop_slab #intensity values in the desired frame
#b_x = bg.x_coords
#b_x_pts = bg.x_loop_pts
#b_y_pts = bg.y_loop_pts
#b_y = bg.y_coords
#b_l = bg.loop_size

# roi coords
#edge = 40
#xmin = int(np.min([f_x_pts[0],f_x_pts[-1]])-edge)
#xmax = int(np.max([f_x_pts[0],f_x_pts[-1]])+edge)
#ymin = int(np.min([f_y_pts[0],f_y_pts[-1]])-edge)
#ymax = int(np.max([f_y_pts[0],f_y_pts[-1]])+edge)

# number of the nodes
if type(nodes_temp)==int: ntemp = nodes_temp
else: ntemp = len(nodes_temp)
if type(nodes_vlos)==int: nvlos = nodes_vlos
else: nvlos = len(nodes_vlos)
if type(nodes_vturb)==int: nvturb = nodes_vturb
else: nvturb = len(nodes_vturb)
if type(nodes_blong)==int: nblong = nodes_blong
else: nblong = len(nodes_blong)
if type(nodes_bhor)==int: nbhor = nodes_bhor
else: nbhor = len(nodes_bhor)
if type(nodes_azi)==int: nazi = nodes_azi
else: nazi = len(nodes_azi)

# smoothing params
fwh = [0.5,1.0,1.5,2.0,2.5]
med = [1,2,3]

# specify depth
dd = 30

for ff in range(len(fwh)):
    for mm in range(len(med)):
        # inversion extract
        f_i = sp.profile(invdir+mode+'_observed_'+index+'.nc') # f_i.dat = (nt, ny, nx, nw, ns)
        f_o = sp.profile(invdir+mode+'_synthetic_cycle'+cyc+'_'+index+'.nc')
        f_m = sp.model(invdir+mode+'_atmosout_cycle'+cyc+'_'+index+'.nc')

        # smoothing
        f_m.smooth(fwhm=fwh[ff], median=med[mm], ntemp=ntemp, nvlos=nvlos, nvturb=nvturb, nBln=nblong, nBho=nbhor, nazi=nazi)
        
        temp = f_m.temp[0,:,:,:].squeeze()*1e-3
        dep = f_m.ltau[0,0,0,:]
        cont = f_i.dat[0,:,:,42,0]
        wave = f_i.wav
        vlos = f_m.vlos[0,:,:,:].squeeze()*1e-5
        vturb = f_m.vturb[0,:,:,:].squeeze()*1e-5
        dy = len(f_m.temp[0,:,0,0])
        dx = len(f_m.temp[0,0,:,0])
        dw = len(wave)
        
        ## figures
        #######

        plt.close('all')
        f = plt.figure(figsize = [6.,5])
        gs = gridspec.GridSpec(2, 1) # grid scale of the 1st col.
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
        ax2.set_xticks(xticks)
        ax2.set_xlabel(r'[pixels]', fontdict=font)
        ax2.set_ylim(0,dy)
        ax2.set_yticks(yticks)
        ax2.set_ylabel(r'[pixels]', fontdict=font)
        
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
        
        # temperature panel
        panel_temp = ax1.imshow(temp[:,:,dd],cmap = 'inferno')
        # overplotting the fibril
        #ax1.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'black',alpha = 0.75,linewidth = 0.5, linestyle = '--')
        #ax1.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'white',alpha = 0.75,linewidth = 0.75, linestyle = '--')
        
        #temp_min, temp_max = np.min(temp[:,:,dd]), np.max(temp[:,:,dd])
        #ct_ticks = np.linspace(np.round(temp_min),np.round(temp_max),np.abs(np.round(temp_min)-np.round(temp_max))+1,dtype='uint')
        ct = plt.colorbar(panel_temp, cax=axin1, orientation="vertical")
        ct.set_label(r'T (kK)', fontdict = font)
        #ct.ax.set_yticklabels(map(str,ct_ticks[:-1])+['>'+str(ct_ticks[-1])], fontdict=font)
        
        # vlos panel
        norm_vlos = np.max(np.abs(vlos[:,:,dd]))
        panel_vlos = ax2.imshow(vlos[:,:,dd], cmap = 'bwr', vmin = -norm_vlos, vmax = norm_vlos)
        #ax2.plot(b_x_pts-xmin,b_y_pts-ymin,color = 'black',alpha = 0.75,linewidth = 0.5, linestyle = '--')
        #ax2.plot(f_x_pts-xmin,f_y_pts-ymin,color = 'black',alpha = 0.75,linewidth = 0.75, linestyle = '--')
        cv = plt.colorbar(panel_vlos, cax=axin2, orientation="vertical")
        cv.set_label(r'v$_{\rm LOS}$ [km/s]', fontdict=font)
        ax1.text(1,dy-21.5,r'FWHM = '+str(fwh[ff])+'\nmedian = '+str(med[mm]),color = 'black',horizontalalignment='left', verticalalignment='bottom', bbox={'facecolor': 'white', 'pad': 1.5,'alpha':0.5})
        plt.show()
        break
        # saving figure
        f.savefig(outdir+'smooth_test/smooth_'+index+'_'+mode+'_fwhm'+str(fwh[ff])+'_median'+str(med[mm])+'.png',dpi = 1000)
        plt.close('all')
        print 'file saved to: ', outdir+'smooth_test/smooth_'+index+'_'+mode+'_fwhm'+str(fwh[ff])+'_median'+str(med[mm])+'.png'
        #break
