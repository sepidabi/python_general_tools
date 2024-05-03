import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


#import spectral as s

#DECLERATIONS
index = '06_165307'
index_bg = '06_165404'

mode = 'roinew_fullstokes'
cyc = '1'
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 9.,
        }

datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
outdir = '/home/seki2695/OUTPUT/inv/'
invdir = '/home/seki2695/INV/stic/roi/'
savedir = datadir+'OUTPUT/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
ss = 0 # stokes param index
res = 0.0375 # CHROMIS pixel size in arcsec
cad = 8.

# CUT coords
x11, x12, x13, x14 = 80, 115, 189, 224
y11, y12, y13, y14 = 6, 42, 110, 147


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
Bln =  f_m.Bln[0,:,:,:].squeeze()*1e-3
dy = len(f_m.temp[0,:,0,0])
dx = len(f_m.temp[0,0,:,0])
dw = len(wave)

# roi coords
edge = 40
xmin, xmax = 1064, 1353 #999+50+20, 1423-20-20-10
ymin, ymax = 487, 639 #457+25+20-30,744-10-70-25 # cropping the ROI

#depths = [52, 28]
depths = [51, 34, 28, 21]

plt.close('all')

#f = plt.figure(figsize = [10.4,7])
f = plt.figure(figsize = [11.5,7.5])
gs = gridspec.GridSpec(4,3) # grid scale of the 1st col.
#plt.style.use('v2.0')
aspect =1# float(f_n_valid)/float(len(f_depRF_clipped))/2.5 #len(f_depRF_clipped)/f_n_valid # panels aspect ratio (y/x)
gs.update(left=0.032,
          right=0.97,
          wspace=0.097,
          bottom=0.06,
          top=0.96,
          hspace = 0.05,
)

# axis info
ytick_pos = np.arange(0,(np.round((ymax-ymin)/40)+1)*40,40)
ytick_lab = ytick_pos*res
xtick_pos = np.arange(0,(np.round((xmax-xmin)/40)+1)*40,40)
xtick_lab = xtick_pos*res

for ii in range(len(depths)):
    dd = depths[ii]

    xticks = np.round(np.linspace(0,280,280/40+1))
    yticks = np.round(np.linspace(0,120,120/30+1))
    ax1 = plt.subplot(gs[ii,0],adjustable = 'box')#,aspect = 'equal')
    ax1.set_xticks([])
    ax1.set_xlim(0,xmax-xmin)
    ax1.set_ylim(0,ymax-ymin)
    ax1.set_yticks(ytick_pos)
    ax1.set_ylabel(r'y [arcsec]', fontdict=font)
    ax1.set_yticklabels(ytick_lab, fontdict = font)
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.tick_params(which='minor', length=2)
    ax1.set_xticks(xtick_pos)
    ax1.set_xticklabels('')

    ax2 = plt.subplot(gs[ii,1],adjustable = 'box')#,aspect = 'equal')
    ax2.set_xticks([])
    ax2.set_xlim(0,xmax-xmin)
    ax2.set_ylim(0,ymax-ymin)
    ax2.set_yticks(ytick_pos)
    ax2.set_yticklabels([])
    ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.tick_params(which='minor', length=2)
    ax2.set_xticks(xtick_pos)
    ax2.set_xticklabels('')
    
    ax3 = plt.subplot(gs[ii,2],adjustable = 'box')#,aspect = 'equal')
    ax3.set_xticks([])
    ax3.set_xlim(0,xmax-xmin)
    ax3.set_ylim(0,ymax-ymin)
    ax3.set_yticks(ytick_pos)
    ax3.set_yticklabels([])
    ax3.set_xlabel('')
    ax3.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax3.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax3.tick_params(which='minor', length=2)
    ax3.set_xticks(xtick_pos)
    ax3.set_xticklabels([])

    if ii ==3:
        ax1.set_xlabel(r'x [arcsec]', fontdict=font)
        ax1.set_xticklabels(xtick_lab, fontdict = font)
        ax2.set_xlabel(r'x [arcsec]', fontdict=font)
        ax2.set_xticklabels(xtick_lab, fontdict = font)
        ax3.set_xlabel(r'x [arcsec]', fontdict=font)
        ax3.set_xticklabels(xtick_lab, fontdict = font)

    else:
        ax1.set_xlabel('')
        ax1.set_xticklabels([])
        ax2.set_xlabel('')
        ax2.set_xticklabels([])
        ax3.set_xlabel('')
        ax3.set_xticklabels([])

    if ii==0:
        ax1.set_title(r'$T$ [kK]', fontdict = font)
        ax2.set_title(r'$v_{\rm LOS}$ [km s$^{-1}$]', fontdict=font)
        ax3.set_title(r'$v_{\rm turb}$ [km s$^{-1}$]', fontdict = font)
        
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
    
        
    if ii==0:
        ct_ticks =np.array( [6.5,7])
        ct_tick_labels = map(str,ct_ticks)
        temp_min, temp_max = np.min(temp[:,:,dd]), ct_ticks[-1]
    elif ii==1:
        ct_ticks =np.array( [4.5, 5.])
        ct_tick_labels = map(str,ct_ticks) #+['>'+str(ct_ticks[-1])]
        temp_min, temp_max = 4.,5.2
    elif ii==2:
        ct_ticks =np.array( [4,5,6])
        ct_tick_labels = map(str,ct_ticks) #+['>'+str(ct_ticks[-1])]
        temp_min, temp_max = 3., 7.   
    else:
        ct_ticks =np.array( [4.5, 5.5, 6.5])
        ct_tick_labels = map(str,ct_ticks) #+['>'+str(ct_ticks[-1])]
        temp_min, temp_max = 4,7
        
    # temperature panel
    #if ii==0:
      #  temp_min, temp_max = np.min(temp[:,:,dd]), np.max(temp[:,:,dd])
    #else:
    #temp_min, temp_max = np.min(temp[:,:,dd]), ct_ticks[-1]    

    print ii, temp_min, temp_max

    # temperature panel
    panel_temp = ax1.imshow(temp[:,:,dd],cmap = 'inferno', vmin = temp_min, vmax = temp_max)

    # magnetic patch contour
    if (ii==0):
        levels =  [-0.8,-0.2,0.2,0.8]
        cp = ax1.contour(Bln[:,:,52], levels = levels, linewidths = 0.95, linestyles = 'solid', alpha = 0.75, colors = ['lime', 'cyan'])
        print ('Bln contour levels', levels)
        #ax1.clabel(cp, inline=True, fontsize=7., fmt = '%1.2f')
    
    ct = plt.colorbar(panel_temp, cax=axin1, orientation="vertical", ticks = ct_ticks)
    ct.ax.set_yticklabels(ct_tick_labels, fontdict = font)
    
    # vlos panel
    sigma_x, sigma_y = 65, 65     # crap in roi
    crap = spnd.filters.gaussian_filter(vlos[:,:,dd], [sigma_x, sigma_y], mode = 'constant')
    vlos_nocrap = vlos[:,:,dd] - crap
    norm_vlos = np.max(np.abs(vlos_nocrap))
    panel_vlos = ax2.imshow(vlos_nocrap, cmap = 'bwr', vmin = -norm_vlos, vmax = norm_vlos)
    cv = plt.colorbar(panel_vlos, cax=axin2, orientation="vertical")
    
    # vturb panel
    norm_vturb = np.max(np.abs(vlos[:,:,dd]))
    panel_vturb = ax3.imshow(vturb[:,:,dd], cmap = 'bone')
    cvmt = plt.colorbar(panel_vturb, cax=axin3, orientation="vertical")

    ax1.text(5,dy-19.5,r'log($\tau_{500}$) = '+str(np.round(dep[dd],decimals = 2)),color = 'black',horizontalalignment='left', verticalalignment='bottom', bbox={'facecolor': 'white', 'pad': 3,'alpha':0.5})

    # fibrils overplot
    need = 1
    if (need):
        lwidth = 0.3 # width of the lines
        alp = 1. #transparency of the plotted lines

        # extracting ALL fibrilS
        for fn in range(len(ffile_fe)):
            # Ca K
            fibdir = datadir+'fr'+str(fr)+'/'
            fib_file = (file_search(fibdir,'crispex*3950*.csav'))[2*fn]#[cat_indx[fn]]
            fib = restore(fibdir+fib_file)
            # Ha
            #fib_file_h = (file_search(fibdir,'crispex*6563*.csav'))[2*fn]#[cat_indx[fn]]
            #fib_h = restore(fibdir+fib_file_h)
            # coords
            f_x = fib.x_coords
            f_x_pts = fib.x_loop_pts
            f_y = fib.y_coords
            f_l = fib.loop_size
            f_y_pts = fib.y_loop_pts

            if  ii==3:
                ax1.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'white', linewidth = lwidth, alpha = alp)
                ax2.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'black', linewidth = lwidth, alpha = alp)
                ax3.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'white', linewidth = lwidth, alpha = alp)

        # extracting THE FIBRIL
        # Ca K
        fibdir = datadir+'fr'+str(fr)+'/'
        fib_file = (file_search(fibdir,'crispex*3950*'+index+'*.csav'))[0]#[cat_indx[fn]]
        fib = restore(fibdir+fib_file)
        # coords
        f_x = fib.x_coords
        f_x_pts = fib.x_loop_pts
        f_y = fib.y_coords
        f_l = fib.loop_size
        f_y_pts = fib.y_loop_pts
        
        if ii==2 or ii==3:
            ax1.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'orangered', linewidth = 1)
            ax2.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'orangered', linewidth = 1)
            ax3.plot(f_x_pts-xmin,f_y_pts-ymin, color = 'orangered', linewidth = 1)



    # cut overplot
    needed_cut = 1
    if (needed_cut):
        if (dd==28):
            # CUT coords
            x11, x12, x13, x14 = 80, 115, 189, 224
            y11, y12, y13, y14 = 6, 42, 110, 147
            
            # point labels
            ax1.text(x14+2, y12-2, 'P', color = 'white', fontsize = 7)
            ax1.text(x13+2, y14-2, 'Q', color = 'white', fontsize = 7)
            ax1.text(x12-7, y11-6, 'R', color = 'white', fontsize = 7)
            ax1.text(x11-7, y13, 'S', color = 'white', fontsize = 7)

            # side cuts
            line_alpha = 0.7
            ax1.plot([x11,x13], [y13,y14], color = 'white', linewidth = 1, alpha = line_alpha)
            ax1.plot([x12,x14], [y11,y12], color = 'white', linewidth = 1, alpha = line_alpha)
            
            # perpen. cuts
            ax1.plot([x12,x11], [y11,y13], color = 'white', linestyle = '--', linewidth = 1, alpha = line_alpha)
            ax1.plot([x14,x13], [y12,y14], color = 'white', linestyle = '--', linewidth = 1, alpha = line_alpha)

plt.show()
filename = outdir+'inv_'+index+'_'+mode+'_cycle'+cyc+'_pathcut_three_new.pdf'
f.savefig(filename, quality = 100)
print 'file saved to: ' + filename
    
