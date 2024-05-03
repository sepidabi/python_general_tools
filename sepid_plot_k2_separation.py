import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
#import spectral as s

fontsize = 10.
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': fontsize,
}

# DECLERATIONS
indxs = ['06_165307', '06_142705','06_165721','07_130930']

plot_mode = 3 # 1 for each fibril profile
                          # 2 for (manual) dI histogram
                          # 3 for seaborn d(I) scatter plot
                          # 4 for seaborn d(T) scatter plot

mode = '_fullstokes'
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

alp = 0.01 #transparency of the plotted lines
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
savedir = datadir+'OUTPUT/'
resdir = '/home/seki2695/INV/stic/'
outdir = '/home/seki2695/OUTPUT/inv/'
pref = ['6302','8542','3950','6563']
pref_name = ['fe','ca8','cak','ha']
fr = 28 #frame no.
ss = 0 # stokes param index
#xx = 100  # set to zero to choose the first pixel
#yy = 0     #stokes I values in the dat array

#maps from the data
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

# calibration info
calib_cak = mf.readfits(datadir+'calib.3934.fits')
calib_ha = mf.readfits(datadir+'calib.6563.fits')

#Mapping the big array of results according to the individual fibril
ffile_fe = file_search(savedir, 'f_obs6302_*.fits')
ffile_ca8 = file_search(savedir, 'f_obs8542_*.fits')
ffile_ck = file_search(savedir, 'f_obs3950_*.fits')
bfile_fe = file_search(savedir, 'b_obs6302_*.fits')
bfile_ca8 = file_search(savedir, 'b_obs8542_*.fits')
bfile_ck = file_search(savedir, 'b_obs3950_*.fits')
#index = ffile_fe
finit = 0
binit = 0
approx_cak_dep_n = 28 # corresponding to log(tau) = -4.59 (RFmax)

# separation containers
f_ksep_blue = np.zeros(0, dtype = np.float64)
f_ksep_red = np.zeros(0, dtype = np.float64)
b_ksep_blue = np.zeros(0, dtype = np.float64)
b_ksep_red = np.zeros(0, dtype = np.float64)

# cak wavelength positions
w_ck = restore(datadir + 'cak_w.sav').spect*1e10 # \AA

for fn in range(len(ffile_fe)):
#for fn in range(10):
    fibril = mf.readfits(savedir + ffile_ck[fn])
    dark = mf.readfits(savedir + bfile_ck[fn])
    fnpix = fibril.shape[2]
    bnpix = dark.shape[2]

    #Fibril extract
    ##========
    amp = 1. #cak amplifier
    #slab
    fibdir = datadir+'fr'+str(fr)+'/'
    fib_file = (file_search(fibdir,'crispex*3950*.csav'))[2*fn]
    fib = restore(fibdir+fib_file)
    f_slab = fib.loop_slab[:,fr,:] #intensity values in the desired frame
    f_x = fib.x_coords
    f_x_pts = fib.x_loop_pts
    f_y = fib.y_coords
    f_l = fib.loop_size
    f_y_pts = fib.y_loop_pts
    #calculating the effective intensity
    #(i.e. the range of wl that contributes to chromospheric brightness)
    f_Itot = np.sum(f_slab[6:15,:])/f_slab[6:15,:].size

    #fibril in Ha
    fib_file_h = (file_search(fibdir,'crispex*6563*.csav'))[2*fn]
    fib_h = restore(fibdir+fib_file_h)
    f_slab_h = fib_h.loop_slab[:,fr,:]/calib_ha[0]  #intensity values in the desired frame
    #calculating the effective intensity
    #(i.e. the range of wl that contributes to chromospheric brightness)
    f_Itot_h = np.sum(f_slab_h[4:10,:])/f_slab_h[4:10,:].size

    #Background extract
    #===========
    #slab
    bg_file = file_search(fibdir,'crispex*3950*.csav')[2*fn+1]
    bg = restore(fibdir+bg_file)
    b_slab = bg.loop_slab[:,fr,:] #intensity values in the desired frame
    b_x = bg.x_coords
    b_x_pts = bg.x_loop_pts
    b_y_pts = bg.y_loop_pts
    b_y = bg.y_coords
    b_l = bg.loop_size
    #calculating the effective intensity
    #(i.e. the range of wl that contributes to chromospheric brightness)
    b_Itot = np.sum(b_slab[6:15,:])/b_slab[6:15,:].size

    # bg in Ha
    bg_file_h = (file_search(fibdir,'crispex*6563*.csav'))[2*fn+1]
    bg_h = restore(fibdir+bg_file_h)
    b_slab_h = bg_h.loop_slab[:,fr,:]/calib_ha[0]  #intensity values in the desired frame
    #calculating the effective intensity
    #(i.e. the range of wl that contributes to chromospheric brightness)
    b_Itot_h = np.sum(b_slab_h[4:10,:])/b_slab_h[4:10,:].size

    # K2 peaks
    c_v, c_r = 9, 10 # k3 rough index
    v1, v2 = 3, 10
    r1,r2 = 9, 16 
    f_k2v = np.max(np.mean(f_slab, axis = 1)[v1:v2]) # k2v fibril
    f_k2v_pos = np.argmax(np.mean(f_slab, axis = 1)[v1:v2]) + v1 # fibril k2v index
    f_k2r = np.max(np.mean(f_slab, axis = 1)[r1:r2]) # k2r fibril
    f_k2r_pos = np.argmax(np.mean(f_slab, axis = 1)[r1:r2]) + r1 # fibril k2r index
    b_k2v = np.max(np.mean(b_slab, axis = 1)[v1:v2]) # k2v bg
    b_k2v_pos = np.argmax(np.mean(b_slab, axis = 1)[v1:v2]) + v1 # k2v bg index
    b_k2r = np.max(np.mean(b_slab, axis = 1)[r1:r2]) # k2r bg
    b_k2r_pos = np.argmax(np.mean(b_slab, axis = 1)[r1:r2]) + r1 # k2r bg index

    #print('fibril no. ', fn)
    # K2 separation
    if (f_k2v > np.mean(f_slab, axis = 1)[c_v] and f_k2r > np.mean(f_slab, axis = 1)[c_r]): # to make sure that the k2 peaks exist
        #print('yay! k2 peak!')
        if (f_k2v > f_k2r):
            #print('blue :D')
            f_ksep_blue = np.append(f_ksep_blue, abs(w_ck[f_k2v_pos] - w_ck[f_k2r_pos]))
            b_ksep_blue = np.append(b_ksep_blue, abs(w_ck[b_k2v_pos] - w_ck[b_k2r_pos]))
        elif (f_k2v < f_k2r):
            #print('red ;)')
            f_ksep_red = np.append(f_ksep_red, abs(w_ck[f_k2v_pos] - w_ck[f_k2r_pos]))
            b_ksep_red = np.append(b_ksep_red, abs(w_ck[b_k2v_pos] - w_ck[b_k2r_pos]))
        else:
            print('equally peaked')
    else:
        print('wierd profile :/','fibril no. ', fn)
    #print('============')


print('STATISTICS:')
print('k2v = ', len(f_ksep_blue)*100/len(ffile_fe), '%')
print('k2r = ', len(f_ksep_red)*100/len(ffile_fe), '%')
print('other profiles ', (len(ffile_fe)-(len(f_ksep_blue)+len(f_ksep_red)))*100/len(ffile_fe), '%')


# whole values
f_ksep = np.concatenate([f_ksep_blue, f_ksep_red])
b_ksep = np.concatenate([b_ksep_blue, b_ksep_red])
    
# the plot
plt.close('all')

f = plt.figure(figsize = (5.5,5))
#ax = plt.subplot(1,1)

xlabel = r' Fibril Average $\rm K_{2}$ Separation $\rm [\AA]$'
ylabel = r' Background Average $\rm K_{2}$ Separation $\rm [\AA]$'

plt.scatter(f_ksep_blue, b_ksep_blue, marker = 's', alpha = 0.04, color = 'b', s = 600, label = 'I$\mathrm{_{K_{v2}}}$ > I$\mathrm{_{K_{r2}}}$')
plt.scatter(f_ksep_red, b_ksep_red, marker = 's', alpha = 0.04, color = 'r', s = 600, label = 'I$\mathrm{_{K_{v2}}}$ < I$\mathrm{_{K_{r2}}}$')
plt.xlabel(xlabel)
plt.ylabel(ylabel)
#plt.text(500,15, 'pearsonr = 0.71; p = 6.7e-29', fontsize = 10.)
plt.xlim(0,1)
plt.ylim(0,1)
plt.legend(markerscale = 0.5)
plt.show()
