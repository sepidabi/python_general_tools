import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import seaborn as sns
import scipy.ndimage as spnd
import matplotlib
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable

# fonts
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }
res = 0.0375 # CHROMIS pixel size in arcsec

# specify the directory and files
resdir = '/home/seki2695/INV/stic/fov/'
dir = resdir + 'results_nocmap/'
outdir = '/home/seki2695/OUTPUT/inv/video/'
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
savedir = datadir+'OUTPUT/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']

# read the fov data from its fits file
header = fits.getheader(resdir + 'fov_inv.fits')
print(header)
fov = fits.getdata(resdir + 'fov_inv.fits')
fov_cmap = np.transpose(fits.getdata(resdir + 'fov_backup_no_cmap.fits'))

# maps
fov_temp = fov[:,:,:,0]
fov_vlos = fov[:,:,:,1]
fov_vlos_cmap = fov_cmap[:,:,:,1]
fov_vturb = fov[:,:,:,2]
fov_Bln = fov[:,:,:,3]
fov_Bln_mod = fov[:,:,:,8]

dep = fov[0,0,:,9]
# dimensions
yy, xx, zz = fov.shape[0], fov.shape[1], fov.shape[2]
fr_i, fr_n = 12, zz-7

# axis info
sc_fact = 266.7 # axis scaling factor
ytick_pos = np.arange(0,(np.round(yy/sc_fact)+1)*sc_fact,sc_fact)
ytick_lab = np.round(ytick_pos*res).astype(int)
xtick_pos = np.arange(0,(np.round(xx/sc_fact)+1)*sc_fact,sc_fact)
xtick_lab = np.round(xtick_pos*res).astype(int)

for dd in range(fr_i, fr_n):

    # Temperature Map
    
    plt.close('all')

    plt.figure(figsize = [1.1*(7.2/2+0.5), 1.1*(8.75-1.75)])
    gs = gridspec.GridSpec(3, 1) # grid scale of the 1st col.

    #plt.style.use('dark_background')
    plt.subplots_adjust(left=0.0,
                        bottom=0.055,
                        right=0.92,
                        top=0.96,
                        wspace=0.2,
                        hspace=0.02
    )
    
    #ax = plt.gca()
    ax_temp = plt.subplot(gs[0,0],adjustable = 'box')
    ax_vlos = plt.subplot(gs[1,0],adjustable = 'box')
    ax_B = plt.subplot(gs[2,0],adjustable = 'box')
    
    temp_panel = unsharp(fov_temp[:,:,dd], alpha = 0.5, sigma = 2.)
    sigma_factor = 3
    temp_min = np.mean(temp_panel) - sigma_factor*np.std(temp_panel)
    temp_max = np.mean(temp_panel) + sigma_factor*np.std(temp_panel)
    im = ax_temp.imshow(temp_panel, cmap = 'gist_heat', origin = 'lower', vmin = temp_min, vmax = temp_max)
    ax_temp.set_title(r'log($\tau _{500}$) = '+str(np.round(dep[dd],decimals = 2)))
    ax_temp.set_xticks([])
    ax_temp.set_xlim(0,xx)
    ax_temp.set_ylim(0,yy)
    ax_temp.set_yticks(ytick_pos)
    ax_temp.set_ylabel(r'y [arcsec]', fontdict=font)
    ax_temp.set_yticklabels(ytick_lab, fontdict = font)
    ax_temp.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax_temp.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax_temp.tick_params(which='minor', length=2)
    ax_temp.set_xticks(xtick_pos)
    ax_temp.set_xticklabels('')
 
    divider = make_axes_locatable(ax_temp)
    #cax = divider.append_axes("right", size="3%", pad=0.1, title = r'$T$ [kK]')
    cax_temp = inset_axes(ax_temp,
                     width="5%",  # width = 5% of parent_bbox width
                     height="100%",  # height : 50%
                     loc='lower left',
                     bbox_to_anchor=(1.02, 0., 1, 1),
                     bbox_transform=ax_temp.transAxes,
                     borderpad=0,
    )

    c_tick = np.round(np.linspace(np.round(temp_min, decimals=0), int(temp_max), int((int(temp_max)-np.round(temp_min, decimals=0))/0.5)+1), decimals = 1)
    cbar = plt.colorbar(im, cax=cax_temp, ticks = c_tick)
    cbar.set_label(r'$T$ [kK]')#, rotation = 0, labelpad = -35, y = 1.05)
    #plt.tight_layout()

    temp_filename = 'fov_temp_dep'
    #plt.savefig(outdir+temp_filename+str(fr_n-(dd)-1).zfill(4)+'.png', dpi = 1000, quality = 100)
    #plt.savefig(outdir+temp_filename+str(fr_n-(dd)-1).zfill(4)+'.pdf', dpi = 1000)


    # vLOS Map
    #plt.style.use('dark_background')
    vlos_panel = unsharp(fov_vlos[:,:,dd], alpha = 0.5, sigma = 2.)
    if dd > 40:
        sigma_factor = 1.5
    else:
        sigma_factor = 2
    vlos_min = -sigma_factor*np.std(np.abs(vlos_panel))
    vlos_max = sigma_factor*np.std(np.abs(vlos_panel))
    im = ax_vlos.imshow(vlos_panel, cmap = 'bwr', origin = 'lower', vmin = vlos_min, vmax = vlos_max)
    ax_vlos.set_xticks([])
    ax_vlos.set_xlim(0,xx)
    ax_vlos.set_ylim(0,yy)
    ax_vlos.set_yticks(ytick_pos)
    ax_vlos.set_yticklabels(ytick_lab, fontdict = font)
    ax_vlos.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax_vlos.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax_vlos.tick_params(which='minor', length=2)
    ax_vlos.set_xticks(xtick_pos)
    ax_vlos.set_xticklabels('')
    ax_vlos.set_ylabel(r'y [arcsec]', fontdict=font)
 
    divider = make_axes_locatable(ax_vlos)
    #cax_vlos = divider.append_axes("right", size="3%", pad=0.1, title = r'$T$ [kK]')
    cax_vlos = inset_axes(ax_vlos,
                     width="5%",  # width = 5% of parent_bbox width
                     height="100%",  # height : 50%
                     loc='lower left',
                     bbox_to_anchor=(1.02, 0., 1, 1),
                     bbox_transform=ax_vlos.transAxes,
                     borderpad=0,
    )
    c_tick = np.round(np.linspace(np.round(vlos_min, decimals=0), int(vlos_max), int((int(vlos_max)-np.round(vlos_min, decimals=0)))+1), decimals = 1)
    cbar = plt.colorbar(im, cax=cax_vlos, ticks = c_tick)
    cbar.set_label(r'$v_{\mathrm{LOS}}$ [km s$^{-1}$] ')#, rotation = 0, labelpad = -35, y = 1.05)
    #plt.tight_layout()

    vlos_filename = 'fov_vlos_dep'
    #plt.savefig(outdir+vlos_filename+str(fr_n-(dd)-1).zfill(4)+'.pdf', dpi = 1000)


    # Bz Map
    #plt.style.use('dark_background')
    Bln_panel = unsharp(fov_Bln[:,:,dd], alpha = 0.5, sigma = 2.)
    sigma_factor = 6
    Bln_min = -sigma_factor*np.std(np.abs(Bln_panel))
    Bln_max = sigma_factor*np.std(np.abs(Bln_panel))
    im = ax_B.imshow(Bln_panel, cmap = 'RdGy_r', origin = 'lower', vmin = Bln_min, vmax = Bln_max)
    ax_B.set_xticks(xtick_pos)
    ax_B.set_xticklabels([])
    ax_B.set_xlim(0,xx)
    ax_B.set_ylim(0,yy)
    ax_B.set_yticks(ytick_pos)
    ax_B.set_yticklabels(ytick_lab, fontdict = font)
    ax_B.set_ylabel(r'y [arcsec]', fontdict=font)
    ax_B.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax_B.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax_B.tick_params(which='minor', length=2)
    ax_B.set_xticks(xtick_pos)
    ax_B.set_xlabel(r'x [arcsec]', fontdict=font)
    ax_B.set_xticklabels(xtick_lab, fontdict = font)
 
    divider = make_axes_locatable(ax_B)
    #cax_B = divider.append_axes("right", size="3%", pad=0.1, title = r'$T$ [kK]')
    cax_B = inset_axes(ax_B,
                     width="5%",  # width = 5% of parent_bbox width
                     height="100%",  # height : 50%
                     loc='lower left',
                     bbox_to_anchor=(1.02, 0., 1, 1),
                     bbox_transform=ax_B.transAxes,
                     borderpad=0,
    )
    c_tick = np.round(np.linspace(np.round(Bln_min, decimals=0), int(Bln_max), int((int(Bln_max)-np.round(Bln_min, decimals=0)))+1), decimals = 1)
    cbar = plt.colorbar(im, cax=cax_B, ticks = c_tick)
    cbar.set_label(r'$B_{\rm long}$ [kG]')#, rotation = 0, labelpad = 35, y = 1.05)

    

    #plt.show()
    #break
    filename = 'fov_all'
    plt.savefig(outdir+filename+str(fr_n-(dd)-1).zfill(4)+'.png', dpi = 1000, quality = 100)

    print('log(tau) = '+str(np.round(dep[dd],decimals = 2))+' -> finished')


# creating the video
os.system("ffmpeg -r 2 -i "+outdir+filename+"%04d.png -vcodec mpeg4 -y "+outdir+filename+".mp4")

# creating the video
#os.system("ffmpeg -r 2 -i "+outdir+vlos_filename+"%04d.png -vcodec mpeg4 -y "+outdir+vlos_filename+".mp4")
