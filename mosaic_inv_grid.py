# This script is meant to mosaic the inversion results
# of the grids of FOV

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


#DECLERATIONS
index = '06_165307'
index_bg = '06_165404'
xmin, xmax, ymin, ymax = 245, 1678, 10, 1130 # FOV coords

# specify the directory and files
resdir = '/home/seki2695/INV/stic/fov/'
dir = resdir + 'results_nocmap/'
outdir = '/home/seki2695/OUTPUT/inv/'
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
savedir = datadir+'OUTPUT/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']

# number of grids, specifeid manually for now
dx, dy = 80, 40 # subregions dimensions
x, y = 17, 28 # number of grids
xx, yy, zz = 1433, 1120, 59 # final mosaic dimensions
fov_temp = np.zeros((zz,xx,yy), dtype = float)
fov_vlos = np.zeros((zz,xx,yy), dtype = float)
fov_vlos_mod = np.zeros((zz,xx,yy), dtype = float)
fov_Bln = np.zeros((zz,xx,yy), dtype = float)
fov_Bln_mod = np.zeros((zz,xx,yy), dtype = float)
fov_Bho = np.zeros((zz,xx,yy), dtype = float)
fov_temp_mod = np.zeros((zz,xx,yy), dtype = float)
fov_vturb = np.zeros((zz,xx,yy), dtype = float)
fov_vturb_mod = np.zeros((zz,xx,yy), dtype = float)
fov_ltau = np.zeros((zz,xx,yy), dtype = float)
dep = 23 # depth to plot upon

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

# other difinitions
res = 0.0375 # CHROMIS pixel size in arcsec
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }

# the grids
for j in range (y):
    for i in range (x):
        grid_dir = 'grid_'+"{0:0=2d}".format(i)+'_'+"{0:0=2d}".format(j)+'/' # grid dir name
        grid = sp.model(dir+grid_dir+'grid_atmosout_cycle1.nc')
        grid_mod = sp.model(dir+grid_dir+'grid_modelin.nc')

        print '===> processing grid_'+"{0:0=2d}".format(i)+'_'+"{0:0=2d}".format(j)
        
        # mosaic dimensions
        x0=i*dx
        x1=(i+1)*dx
        y0=j*dy
        y1=(j+1)*dy

        # maps
        fov_temp[:, x0:x1, y0:y1] = np.transpose(grid.temp[0,:,:,:])
        fov_temp_mod[:, x0:x1, y0:y1] = np.transpose(grid_mod.temp[0,:,:,:])
        fov_vlos[:, x0:x1, y0:y1] = np.transpose(grid.vlos[0,:,:,:])
        fov_vlos_mod[:, x0:x1, y0:y1] = np.transpose(grid_mod.vlos[0,:,:,:])
        fov_vturb[:, x0:x1, y0:y1] = np.transpose(grid.vturb[0,:,:,:])
        fov_vturb_mod[:, x0:x1, y0:y1] = np.transpose(grid_mod.vturb[0,:,:,:])
        fov_Bln[:, x0:x1, y0:y1] = np.transpose(grid.Bln[0,:,:,:])
        fov_Bln_mod[:, x0:x1, y0:y1] = np.transpose(grid_mod.Bln[0,:,:,:])
        fov_Bho[:, x0:x1, y0:y1] = np.transpose(grid.Bho[0,:,:,:])
        fov_ltau[:, x0:x1, y0:y1] = np.transpose(grid.ltau[0,:,:,:])
        tau = grid.ltau[0,0,0,:]

        
    grid_ver_dir  = 'grid_ver_'+"{0:0=2d}".format(j)+'/' # grid dir name
    grid_ver = sp.model(dir+grid_ver_dir+'grid_atmosout_cycle1.nc')
    grid_ver_mod = sp.model(dir+grid_ver_dir+'grid_modelin.nc')
    print '===> processing grid_ver_'+"{0:0=2d}".format(j)
    
    # filling the maps
    fov_temp[:, x*dx:, y0:y1] = np.transpose(grid_ver.temp[0,:,:,:])
    fov_temp_mod[:, x*dx:, y0:y1] = np.transpose(grid_ver_mod.temp[0,:,:,:])
    fov_vlos[:, x*dx:, y0:y1] = np.transpose(grid_ver.vlos[0,:,:,:])
    fov_vlos_mod[:, x*dx:, y0:y1] = np.transpose(grid_ver_mod.vlos[0,:,:,:])
    fov_vturb[:, x*dx:, y0:y1] = np.transpose(grid_ver.vturb[0,:,:,:])
    fov_vturb_mod[:, x*dx:, y0:y1] = np.transpose(grid_ver_mod.vturb[0,:,:,:])
    fov_Bln[:, x*dx:, y0:y1] = np.transpose(grid_ver.Bln[0,:,:,:])
    fov_Bln_mod[:, x*dx:, y0:y1] = np.transpose(grid_ver_mod.Bln[0,:,:,:])
    fov_Bho[:, x*dx:, y0:y1] = np.transpose(grid_ver.Bho[0,:,:,:])
    
# fov cube container
fov = np.zeros((10,59,1433,1120), dtype=float)

# filling the fov values
fov[0,:,:] = fov_temp.squeeze()*1e-3
fov[5,:,:] = fov_temp_mod.squeeze()*1e-3
fov[1,:,:] = fov_vlos.squeeze()*1e-5
fov[6,:,:] = fov_vlos_mod.squeeze()*1e-5
fov[2,:,:] = fov_vturb.squeeze()*1e-5
fov[7,:,:] = fov_vturb_mod.squeeze()*1e-5
fov[3,:,:] = fov_Bln.squeeze()*1e-3
fov[8,:,:] = fov_Bln_mod.squeeze()*1e-3
fov[4,:,:] = fov_Bho.squeeze()*1e-3
fov[9,:,:] = fov_ltau.squeeze()*1.

# transposing to make the dimensions compatible with observations
fov = np.transpose(fov)
fov_no_cmap = fov

# apply cmap correction
if(1):
    for tt in range(fov.shape[2]):
        sigma_x, sigma_y = 65, 65
        crap = spnd.filters.gaussian_filter(fov[:,:,tt,1], [sigma_x, sigma_y], mode = 'constant')
        crap_NN = spnd.filters.gaussian_filter(fov[:,:,tt,6], [sigma_x, sigma_y], mode = 'constant')
        fov[:,:,tt,1] = fov[:,:,tt,1] - crap
        fov[:,:,tt,6] = fov[:,:,tt,6] - crap_NN
        print(tau[tt], crap)

# save the fits file
if (0):
    # creating the fits file of the fov inv results
    hdu = fits.PrimaryHDU(data = fov)
    hdu_no_cmap = fits.PrimaryHDU(data = fov_no_cmap)
    
    # modify the header
    hdu.header.set('NAXIS', comment = 'shape = ' + str(fov.shape))
    hdu.header.rename_keyword('NAXIS1', 'y-axis')
    hdu.header.set('y-axis', comment = 'pxls')
    hdu.header.rename_keyword('NAXIS2', 'x-axis')
    hdu.header.set('x-axis', comment = 'pxls')
    hdu.header.rename_keyword('NAXIS3', 'logt')
    hdu.header.rename_keyword('NAXIS4', 'params')
    
    # 'T[kK], vLOS[km/s], vturb[km/s], Bln[kG], Bho[kG], T_NN, vLOS_NN, vturb_NN, Bln_NN, tau'
    hdu.header['param0'] = 'T [kK]'
    hdu.header['param1'] = 'vLOS[km/s]'
    hdu.header['param2'] = 'vturb[km/s]'
    hdu.header['param3'] = 'Bln[kG]'
    hdu.header['param4'] = 'Bho[kG]'
    hdu.header['param5'] = 'T_NN'
    hdu.header['param6'] = 'vLOS_NN'
    hdu.header['param7'] = 'vturb_NN'
    hdu.header['param8'] = 'Bln_NN'
    hdu.header['param9'] = 'tau'

    # saving...
    hdu.writeto(resdir + 'fov_inv.fits', overwrite = True)
    hdu_no_cmap.writeto(resdir + 'fov_inv_no_cmap.fits', overwrite = True)
    
