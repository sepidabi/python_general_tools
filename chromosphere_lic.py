# A brief script to generate cool images
# for the thesis
# from chromospheric observation at SST
# written by Sepideh Kianfar

from sepid import *

datadir = '/scratch/jaime/CHROMIS/2016.09.19_fluxem_09_28_36_cakfine_crisp/'
outdir = '/home/seki2695/OUTPUT/'
w_cak = 11
w_ha = 7
fr = 16

cak_cube = datadir + file_search(datadir, 'crispex_3950*time-corrected.fcube')[0]
cak_cube_sp = datadir + file_search(datadir, 'crispex_3950*time-corrected_sp.fcube')[0]
#nx, ny, dum, ns, dtype, ndim = lp.lphead(cak_cube)
#nw_cak, nt, dum, ns, dtype, ndim = lp.lphead(cak_cube_sp)

ha_cube = datadir + file_search(datadir, 'crispex.6563*time_corrected_CHROMIS.fcube')[0]
ha_cube_sp = datadir + file_search(datadir, 'crispex.6563*time_corrected_CHROMIS_sp.fcube')[0]
#nx, ny, dum, ns, dtype, ndim = lp.lphead(cak_cube)
#nw_ha, nt, dum, ns, dtype, ndim = lp.lphead(cak_cube_sp)

cak_scan = unsharp(lp_read(cak_cube, cak_cube_sp)[fr, 0, w_cak, :, :], alpha = 0.5, sigma = 2.)

ha_scan = unsharp(lp_read(ha_cube, ha_cube_sp)[fr, 0, w_ha, :, :], alpha = 0.5, sigma = 2.)

# Ca II K intensity
plt.close('all')

# clipping the image
sigma_factor = 4.
cak_min = np.mean(cak_scan) - 2.5*np.std(cak_scan)
cak_max = np.mean(cak_scan) + 6*np.std(cak_scan)

# plotting the map
plt.imshow(cak_scan, origin = 'lower', cmap = 'gray', vmin = cak_min, vmax = cak_max)
plt.xticks([])
plt.yticks([])

plt.tight_layout()
plt.show()
filename = 'chromosphere_cak.pdf'
plt.savefig(outdir + filename, dpi = 1000)

# Ha intensity
plt.close('all')

# clipping the image
sigma_factor = 4.
ha_min = np.mean(ha_scan) - 2.5*np.std(ha_scan)
ha_max = np.mean(ha_scan) + 6*np.std(ha_scan)

# plotting the map
plt.imshow(ha_scan, origin = 'lower', cmap = 'gray', vmin = ha_min, vmax = ha_max)
plt.xticks([])
plt.yticks([])

plt.tight_layout()
plt.show()
filename = 'chromosphere_ha.pdf'
plt.savefig(outdir + filename, dpi = 1000)
