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
        'color':  'white',
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
    plt.style.use('dark_background')
    plt.figure(figsize = [13*0.6,10*0.6])
    plt.subplots_adjust(left=0.07,
                        bottom=0.,
                        right=0.88,
                        top=1.,
                        wspace=0.2,
                        hspace=0.2
    )
    
    ax = plt.gca()
    temp_panel = unsharp(fov_temp[:,:,dd], alpha = 0.5, sigma = 2.)
    sigma_factor = 3
    temp_min = np.mean(temp_panel) - sigma_factor*np.std(temp_panel)
    temp_max = np.mean(temp_panel) + sigma_factor*np.std(temp_panel)
    im = plt.imshow(temp_panel, cmap = 'gist_heat', origin = 'lower', vmin = temp_min, vmax = temp_max)
    ax.set_title(r'log($\tau _{500}$) = '+str(np.round(dep[dd],decimals = 2)))
    ax.set_xlim(0,xx)
    ax.set_ylim(0,yy)
    ax.set_yticks(ytick_pos)
    ax.set_ylabel(r'y [arcsec]', fontdict=font)
    ax.set_yticklabels(ytick_lab, fontdict = font)
    ax.set_xticks(xtick_pos)
    ax.set_xlabel(r'x [arcsec]', fontdict=font)
    ax.set_xticklabels(xtick_lab, fontdict = font)
 
    divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="3%", pad=0.1, title = r'$T$ [kK]')
    cax = inset_axes(ax,
                     width="5%",  # width = 5% of parent_bbox width
                     height="100%",  # height : 50%
                     loc='lower left',
                     bbox_to_anchor=(1.02, 0., 1, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0,
    )

    c_tick = np.round(np.linspace(np.round(temp_min, decimals=0), int(temp_max), int((int(temp_max)-np.round(temp_min, decimals=0))/0.5)+1), decimals = 1)
    cbar = plt.colorbar(im, cax=cax, ticks = c_tick)
    cbar.set_label(r'$T$ [kK]', rotation = 0, labelpad = -35, y = 1.05)
    #plt.tight_layout()

    temp_filename = 'fov_temp_dep'
    plt.savefig(outdir+temp_filename+str(fr_n-(dd)-1).zfill(4)+'.png', dpi = 1000, quality = 100)
    plt.savefig(outdir+temp_filename+str(fr_n-(dd)-1).zfill(4)+'.pdf', dpi = 1000)


    # vLOS Map
    plt.close('all')
    plt.style.use('dark_background')
    plt.figure(figsize = [13*0.6,10*0.6])
    plt.subplots_adjust(left=0.07,
                        bottom=0.,
                        right=0.87,
                        top=1.,
                        wspace=0.2,
                        hspace=0.2
    )
    
    ax = plt.gca()
    vlos_panel = unsharp(fov_vlos[:,:,dd], alpha = 0.5, sigma = 2.)
    if dd > 40:
        sigma_factor = 1.5
    else:
        sigma_factor = 2
    vlos_min = -sigma_factor*np.std(np.abs(vlos_panel))
    vlos_max = sigma_factor*np.std(np.abs(vlos_panel))
    im = plt.imshow(vlos_panel, cmap = 'bwr', origin = 'lower', vmin = vlos_min, vmax = vlos_max)
    ax.set_title(r'log($\tau _{500}$) = '+str(np.round(dep[dd],decimals = 2)))
    ax.set_xlim(0,xx)
    ax.set_ylim(0,yy)
    ax.set_yticks(ytick_pos)
    ax.set_ylabel(r'y [arcsec]', fontdict=font)
    ax.set_yticklabels(ytick_lab, fontdict = font)
    ax.set_xticks(xtick_pos)
    ax.set_xlabel(r'x [arcsec]', fontdict=font)
    ax.set_xticklabels(xtick_lab, fontdict = font)
 
    divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="3%", pad=0.1, title = r'$T$ [kK]')
    cax = inset_axes(ax,
                     width="5%",  # width = 5% of parent_bbox width
                     height="100%",  # height : 50%
                     loc='lower left',
                     bbox_to_anchor=(1.02, 0., 1, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0,
    )
    c_tick = np.round(np.linspace(np.round(vlos_min, decimals=0), int(vlos_max), int((int(vlos_max)-np.round(vlos_min, decimals=0)))+1), decimals = 1)
    cbar = plt.colorbar(im, cax=cax, ticks = c_tick)
    cbar.set_label(r'$v_{\mathrm{LOS}}$ [km s$^{-1}$] ', rotation = 0, labelpad = -35, y = 1.05)
    #plt.tight_layout()

    vlos_filename = 'fov_vlos_dep'
    plt.savefig(outdir+vlos_filename+str(fr_n-(dd)-1).zfill(4)+'.png', dpi = 1000, quality = 100)
    plt.savefig(outdir+vlos_filename+str(fr_n-(dd)-1).zfill(4)+'.pdf', dpi = 1000)
    
    print('log(tau) = '+str(np.round(dep[dd],decimals = 2))+' -> finished')


# creating the video
os.system("ffmpeg -r 2 -i "+outdir+temp_filename+"%04d.png -vcodec mpeg4 -y "+outdir+temp_filename+".mp4")

# creating the video
os.system("ffmpeg -r 2 -i "+outdir+vlos_filename+"%04d.png -vcodec mpeg4 -y "+outdir+vlos_filename+".mp4")
