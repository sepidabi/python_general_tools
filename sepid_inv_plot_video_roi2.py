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

mode = 'roinew_fullstokes'
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
savedir = datadir+'OUTPUT/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
ss = 0 # stokes param index
res = 0.0375 # CHROMIS pixel size in arcsec


# data extract
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
fr = 28 #frame no.
file_fe =file_search(datadir,'crispex*'+pref[0]+'*.fcube')
file_ca8 = file_search(datadir,'crispex*'+pref[1]+'*.fcube')
file_cak = file_search(datadir,'crispex*'+pref[2]+'*.fcube')
file_ha = file_search(datadir,'crispex*'+pref[3]+'*.fcube')
cube_fe = lp_read(datadir+file_fe[0],datadir+file_fe[1])[fr]
cube_ca8 = lp_read(datadir+file_ca8[0],datadir+file_ca8[1])[fr]
cube_cak = lp_read(datadir+file_cak[0],datadir+file_cak[1])[fr]
cube_ha = lp_read(datadir+file_ha[0],datadir+file_ha[1])[fr]
c_map = cube_cak[0,-1,:,:]
cak_int_un = (lp_read_scan(datadir+'cak_int_un.fcube'))[fr,:,:]

# fibril calibrated files
ffile_fe = file_search(savedir, 'f_obs6302_*.fits')
bfile_fe = file_search(savedir, 'b_obs6302_*.fits')

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

# roi coords
edge = 40
xmin, xmax = 1064, 1353 #999+50+20, 1423-20-20-10
ymin, ymax = 487, 639 #457+25+20-30,744-10-70-25 # cropping the ROI

# preparing for the video
#fr_i, fr_n = 5, len(dep)-5
fr_i, fr_n = 12, len(dep)-7

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

    # axis info
    ytick_pos = np.arange(0,(np.round((ymax-ymin)/40)+1)*40,40)
    ytick_lab = ytick_pos*res
    xtick_pos = np.arange(0,(np.round((xmax-xmin)/40)+1)*40,40)
    xtick_lab = xtick_pos*res

    #xticks = np.round(np.linspace(0,280,280/40+1))
    #yticks = np.round(np.linspace(0,120,120/30+1))
    ax1 = plt.subplot(gs[0,0],adjustable = 'box')#,aspect = 'equal')
    ax1.set_xticks(xtick_pos)
    ax1.set_xticklabels('')
    ax1.set_xlim(0,dx)
    ax1.set_ylim(0,dy)
    ax1.set_yticks(ytick_pos)
    ax1.set_ylabel(r'y [arcsec]', fontdict=font)
    ax1.set_yticklabels(ytick_lab)
    ax2 = plt.subplot(gs[1,0],adjustable = 'box')#,aspect = 'equal')
    ax2.set_xticks(xtick_pos)
    ax2.set_xticklabels('')
    ax2.set_xlim(0,dx)
    ax2.set_ylim(0,dy)
    ax2.set_yticks(ytick_pos)
    ax2.set_yticklabels(ytick_lab)
    ax2.set_ylabel(r'y [arcsec]', fontdict=font)
    ax3 = plt.subplot(gs[2,0],adjustable = 'box')#,aspect = 'equal')
    ax3.set_xlabel(r'x [arcsec]', fontdict=font)
    ax3.set_xlim(0,dx)
    ax3.set_ylim(0,dy)
    ax3.set_yticks(ytick_pos)
    ax3.set_yticklabels(ytick_lab)
    ax3.set_xticks(xtick_pos)
    ax3.set_xticklabels(xtick_lab)
    ax3.set_ylabel(r'y [arcsec]', fontdict=font)
    
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
    ct = plt.colorbar(panel_temp, cax=axin1, orientation="vertical")
    ct.set_label(r'$T$ [kK]', fontdict = font)
    
    # vlos panel
    sigma_x, sigma_y = 65, 65     # crap in roi
    crap = spnd.filters.gaussian_filter(vlos[:,:,dd], [sigma_x, sigma_y], mode = 'constant')
    vlos_nocrap = vlos[:,:,dd] - crap
    norm_vlos = np.max(np.abs(vlos_nocrap))-1
    panel_vlos = ax2.imshow(vlos_nocrap, cmap = 'bwr', vmin = -norm_vlos, vmax = norm_vlos)
    cv = plt.colorbar(panel_vlos, cax=axin2, orientation="vertical")
    cv.set_label(r'$v_{\rm LOS}$ [km s$^{-1}$]', fontdict=font)
    
    # vturb panel
    panel_vturb = ax3.imshow(vturb[:,:,dd], cmap = 'bone')
    cvmt = plt.colorbar(panel_vturb, cax=axin3, orientation="vertical")
    cvmt.set_label(r'$v_{\rm turb}$ [km s$^{-1}$]', fontdict = font)
    ax1.text(5,dy-17,r'log($\tau_{500}$) = '+str(np.round(dep[dd],decimals = 2)),color = 'black',horizontalalignment='left', verticalalignment='bottom', bbox={'facecolor': 'white', 'pad': 3,'alpha':0.75})

    '''
    if (dd==23 or dd==28):

       for fn in range(len(ffile_fe)):
            # extracting FIBRIL
            # Ca K
            fibdir = datadir+'fr'+str(fr)+'/'
            fib_file = (file_search(fibdir,'crispex*3950*.csav'))[2*fn]#[cat_indx[fn]]
            fib = restore(fibdir+fib_file)
            # Ha
            fib_file_h = (file_search(fibdir,'crispex*6563*.csav'))[2*fn]#[cat_indx[fn]]
            fib_h = restore(fibdir+fib_file_h)
            # coords
            f_x = fib.x_coords
            f_x_pts = fib.x_loop_pts
            f_y = fib.y_coords
            f_l = fib.loop_size
            f_y_pts = fib.y_loop_pts
            
            # extracting BG
            bg_file = file_search(fibdir,'crispex*3950*.csav')[2*fn+1]#[cat_indx[fn]+1]
            bg = restore(fibdir+bg_file)
            # Ha
            bg_file_h = (file_search(fibdir,'crispex*6563*.csav'))[2*fn+1]#[cat_indx[fn]+1]
            bg_h = restore(fibdir+bg_file_h)
            
            # coords
            b_x = bg.x_coords
            b_x_pts = bg.x_loop_pts
            b_y_pts = bg.y_loop_pts
            b_y = bg.y_coords
            b_l = bg.loop_size
            
            lwidth = 1.5 # width of the lines
            alp = 0.5 #transparency of the plotted lines
            
            ax1.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'red', linewidth = lwidth, alpha = alp, linestyle = '-')
            ax1.plot(b_x_pts-xmin,b_y_pts-ymin, color = 'grey', linewidth = lwidth, alpha = alp, linestyle = '-')
            ax2.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'red', linewidth = lwidth, alpha = alp, linestyle = '-')
'''
    # to make the video starting from the photosphere
    f.savefig(outdir+'inv_'+index+'_'+mode+'_cycle'+cyc+'dep'+str(fr_n-(dd)-1).zfill(4)+'.png', dpi = 1000)
    

'''    if (dd==28):
        # CUT coords
        x11, x12, x13, x14 = 80, 115, 189, 224
        y11, y12, y13, y14 = 6, 42, 110, 147
        
        # side cuts
        ax1.plot([x11,x13], [y13,y14], color = 'white', linewidth = 1)
        ax1.plot([x12,x14], [y11,y12], color = 'white', linewidth = 1)
        
        # perpen. cuts
        ax1.plot([x12,x11], [y11,y13], color = 'white', linestyle = '--', linewidth = 1)
        ax1.plot([x14,x13], [y12,y14], color = 'white', linestyle = '--', linewidth = 1)
            
        f.savefig(outdir+'inv_'+index+'_'+mode+'_cycle'+cyc+'dep'+str(fr_n-(dd)-1).zfill(4)+'_pathcut.pdf', quality = 100)
    
        plt.show()
'''            #break

# creating the video
os.system("ffmpeg -r 2 -i "+outdir+'inv_'+index+'_'+mode+'_cycle'+cyc+'dep'+"%04d.png -vcodec mpeg4 -y "+outdir+'video/inv_'+index+'_'+mode+'_fullstokes_cycle'+cyc+"_new.mp4")
