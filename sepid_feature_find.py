import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
from matplotlib.gridspec import GridSpec

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

# k2 separation containers
f_ksep_blue = np.zeros(0, dtype = np.float64)
f_ksep_red = np.zeros(0, dtype = np.float64)
b_ksep_blue = np.zeros(0, dtype = np.float64)
b_ksep_red = np.zeros(0, dtype = np.float64)
f_ksep = np.zeros(0, dtype = np.float64)
b_ksep = np.zeros(0, dtype = np.float64)
f_A = np.zeros(0, dtype = np.float64)
b_A = np.zeros(0, dtype = np.float64)
fI_k2v = np.zeros(0, dtype = np.float64)
fI_k2r = np.zeros(0, dtype = np.float64)
bI_k2v = np.zeros(0, dtype = np.float64)
bI_k2r = np.zeros(0, dtype = np.float64)
px_ntot = 0.
f_ntot = 0.
f_b_ntot = 0.

# k1 separation containers
f_k1sep_blue = np.zeros(0, dtype = np.float64)
f_k1sep_red = np.zeros(0, dtype = np.float64)
b_k1sep_blue = np.zeros(0, dtype = np.float64)
b_k1sep_red = np.zeros(0, dtype = np.float64)
f_k1sep = np.zeros(0, dtype = np.float64)
b_k1sep = np.zeros(0, dtype = np.float64)
f_Ak1 = np.zeros(0, dtype = np.float64)
b_Ak1 = np.zeros(0, dtype = np.float64)
fI_k1v = np.zeros(0, dtype = np.float64)
fI_k1r = np.zeros(0, dtype = np.float64)
bI_k1v = np.zeros(0, dtype = np.float64)
bI_k1r = np.zeros(0, dtype = np.float64)
fI_k1k2v = np.zeros(0, dtype = np.float64)
fI_k1k2r = np.zeros(0, dtype = np.float64)
bI_k1k2v = np.zeros(0, dtype = np.float64)
bI_k1k2r = np.zeros(0, dtype = np.float64)
f_k1ntot = 0.
f_b_k1ntot = 0.


pixels = 0 # total number of pixels with mustache profile

# cak wavelength positions
w_ck = restore(datadir + 'cak_w.sav').spect*1.e10 # \AA

plt.close('all')

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
    f_slab = fib.loop_slab[:,fr,:]*1000  #intensity values in the desired frame
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
    b_slab = bg.loop_slab[:,fr,:]*1000 #intensity values in the desired frame
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
    v1, v2 = 3, 9
    r1,r2 = 10, 17
    for pp in range(np.min([f_slab_h.shape[1], b_slab_h.shape[1]])):

        # total number of fibrillar pixels
        px_ntot+=1

        # K2
        f_k2v = np.max(f_slab[v1:v2,pp]) # k2v fibril
        f_k2v_pos = np.argmax(f_slab[v1:v2,pp]) + v1 # fibril k2v index
        f_k2r = np.max(f_slab[r1:r2,pp]) # k2r fibril
        f_k2r_pos = np.argmax(f_slab[r1:r2,pp]) + r1 # fibril k2r index
        b_k2v = np.max(b_slab[v1:v2,pp]) # k2v bg
        b_k2v_pos = np.argmax(b_slab[v1:v2,pp]) + v1 # k2v bg index
        b_k2r = np.max(b_slab[r1:r2,pp]) # k2r bg
        b_k2r_pos = np.argmax(b_slab[r1:r2,pp]) + r1 # k2r bg index

        # K1
        f_k1v = np.min(f_slab[0:f_k2v_pos,pp]) # k1v fibril
        f_k1v_pos = np.argmin(f_slab[0:f_k2v_pos,pp]) # fibril k1v index
        f_k1r = np.min(f_slab[f_k2r_pos:20,pp]) # k1r fibril
        f_k1r_pos = np.argmin(f_slab[f_k2r_pos:20,pp]) + f_k2r_pos # fibril k1r index
        b_k1v = np.min(b_slab[0:f_k2v_pos,pp]) # k1v bg
        b_k1v_pos = np.argmin(b_slab[0:f_k2v_pos,pp]) # k1v bg index
        b_k1r = np.min(b_slab[f_k2r_pos:20,pp]) # k1r bg
        b_k1r_pos = np.argmin(b_slab[f_k2r_pos:20,pp]) + f_k2r_pos # k1r bg index

        # double-peaked fibrils
        if (f_k2v > f_slab[f_k2v_pos-1, pp] and
            f_k2v > f_slab[f_k2v_pos+1, pp] and
            f_k2r > f_slab[f_k2r_pos-1, pp] and
            f_k2r > f_slab[f_k2r_pos+1, pp] and
            f_k1v < f_slab[f_k1v_pos-1, pp] and
            f_k1v < f_slab[f_k1v_pos+1, pp] and
            f_k1r < f_slab[f_k1r_pos-1, pp] and
            f_k1r < f_slab[f_k1r_pos+1, pp]
        ):


            f_ntot+=1
            # double-peaked in both fb and bg
            if(b_k2v > b_slab[b_k2v_pos-1, pp] and
               b_k2v > b_slab[b_k2v_pos+1, pp] and
               b_k2r > b_slab[b_k2r_pos-1, pp] and
               b_k2r > b_slab[b_k2r_pos+1, pp] and
               b_k1v < b_slab[b_k1v_pos-1, pp] and
               b_k1v < b_slab[b_k1v_pos+1, pp] and
               b_k1r < b_slab[b_k1r_pos-1, pp] and
               b_k1r < b_slab[b_k1r_pos+1, pp]
            ):

                f_b_ntot+=1


                # K2 Asymmetry
                fI_k2v, fI_k2r = np.append(fI_k2v, f_k2v), np.append(fI_k2r, f_k2r)
                bI_k2v, bI_k2r = np.append(bI_k2v, b_k2v), np.append(bI_k2r, b_k2r)
                f_A = np.append(f_A, (f_k2v - f_k2r)/(f_k2v + f_k2r))
                b_A = np.append(b_A, (b_k2v - b_k2r)/(b_k2v + b_k2r))
                
                # K2 separation
                f_ksep = np.append(f_ksep, abs(w_ck[f_k2v_pos] - w_ck[f_k2r_pos]))
                b_ksep = np.append(b_ksep, abs(w_ck[b_k2v_pos] - w_ck[b_k2r_pos]))                        
                
                if (f_k2v > f_k2r):
                    #print('blue :D')
                    f_ksep_blue = np.append(f_ksep_blue, abs(w_ck[f_k2v_pos] - w_ck[f_k2r_pos]))
                    b_ksep_blue = np.append(b_ksep_blue, abs(w_ck[b_k2v_pos] - w_ck[b_k2r_pos]))
                elif (f_k2v < f_k2r):
                    #print('red ;)')
                    f_ksep_red = np.append(f_ksep_red, abs(w_ck[f_k2v_pos] - w_ck[f_k2r_pos]))
                    b_ksep_red = np.append(b_ksep_red, abs(w_ck[b_k2v_pos] - w_ck[b_k2r_pos]))

                                    # K1 Asymmetry
                fI_k1v, fI_k1r = np.append(fI_k1v, f_k1v), np.append(fI_k1r, f_k1r)
                bI_k1v, bI_k1r = np.append(bI_k1v, b_k1v), np.append(bI_k1r, b_k1r)
                f_Ak1 = np.append(f_Ak1, (f_k1v - f_k1r)/(f_k1v + f_k1r))
                b_Ak1 = np.append(b_Ak1, (b_k1v - b_k1r)/(b_k1v + b_k1r))
                
                # K1 separation
                f_k1sep = np.append(f_k1sep, abs(w_ck[f_k1v_pos] - w_ck[f_k1r_pos]))
                b_k1sep = np.append(b_k1sep, abs(w_ck[b_k1v_pos] - w_ck[b_k1r_pos]))
                
                #print('yay! k1 peak!')
                if (f_k1v > f_k1r):
                    #print('blue :D')
                    f_k1sep_blue = np.append(f_k1sep_blue, abs(w_ck[f_k1v_pos] - w_ck[f_k1r_pos]))
                    b_k1sep_blue = np.append(b_k1sep_blue, abs(w_ck[b_k1v_pos] - w_ck[b_k1r_pos]))
                elif (f_k1v < f_k1r):
                    #print('red ;)')
                    f_k1sep_red = np.append(f_k1sep_red, abs(w_ck[f_k1v_pos] - w_ck[f_k1r_pos]))
                    b_k1sep_red = np.append(b_k1sep_red, abs(w_ck[b_k1v_pos] - w_ck[b_k1r_pos]))
                else:
                    print('equally peaked')


            else:
                print('wierd profile :/','fibril no. ', fn, ' - pixel ', pp)
                #print('============')
                
    pixels +=np.min([f_slab_h.shape[1], b_slab_h.shape[1]])

# whole values
f_ksep = np.concatenate([f_ksep_blue, f_ksep_red])
b_ksep = np.concatenate([b_ksep_blue, b_ksep_red])

# percentages
print('STATISTICS:')
k2v = 100.*(np.where(f_A>0))[0].size/f_b_ntot #f_A.size
k2r = 100.*(np.where(f_A<0))[0].size/f_b_ntot #f_A.size
#k2 = 100. - k2v - k2r
print(f_ntot*100/px_ntot, '% 2-peak fibrillar pixels:')
print(f_b_ntot*100/f_ntot, '% 2-peak fb and bg pixels')
print(k2v, '% K2V')
print(k2r, '% K2R')
print(100-k2v-k2r, 'equally peaked')



################
# the K2-ASYMETRY plot
################
plt.close('all')
fig = plt.figure(figsize = (5,5))

# k1 and k2 intensitty differences

x = 100*(np.mean([fI_k2v,fI_k2r], axis = 0, dtype=float) - np.mean([bI_k2v,bI_k2r], axis = 0, dtype=float))/np.mean([bI_k2v,bI_k2r], axis = 0, dtype=float)

y = 100*(np.mean([fI_k1v,fI_k1r], axis = 0, dtype=float) - np.mean([bI_k1v,bI_k1r], axis = 0, dtype=float))/np.mean([bI_k1v,bI_k1r], axis = 0, dtype=float)

xmin, xmax = -15,35
ymin, ymax = -25,30

gs = GridSpec(5,5)

ax_joint = fig.add_subplot(gs[1:5,0:4])
ax_marg_x = fig.add_subplot(gs[0,0:4])
ax_marg_y = fig.add_subplot(gs[1:5,4])

# Turn off ticks and labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
plt.setp(ax_marg_y.get_xticklabels(), visible=False)
plt.setp(ax_marg_x.get_yticklabels(), visible=False)
ax_marg_x.get_yaxis().set_visible(False)
ax_marg_y.get_xaxis().set_visible(False)


# horizontal hist
ax_marg_x.spines['top'].set_visible(False)
ax_marg_x.spines['right'].set_visible(False)
ax_marg_x.spines['left'].set_visible(False)
# vertical hist
ax_marg_y.spines['top'].set_visible(False)
ax_marg_y.spines['right'].set_visible(False)
ax_marg_y.spines['bottom'].set_visible(False)
# joint plot
ax_joint.spines['top'].set_visible(False)
ax_joint.spines['right'].set_visible(False)

# Set labels on joint
ax_joint.set_xlabel(r'$\mathrm{\Delta I _{K_{2}}}$ [%]', fontdict = font)
ax_joint.set_ylabel(r'$\mathrm{\Delta I _{K_{1}}}$ [%]', fontdict = font)
sns.kdeplot(x,y, cmap = 'Reds', ax = ax_joint, shade = True)
sns.kdeplot(x, ax = ax_marg_x, shade=True, color = 'orangered')
sns.kdeplot(y, ax = ax_marg_y, shade=True, color = 'orangered', vertical = True)

# Set ax limits on marginals
ax_joint.set_xlim([xmin,xmax])
ax_joint.set_ylim([ymin,ymax])
ax_marg_y.set_ylim([ymin,ymax])
ax_marg_x.set_xlim([xmin,xmax])
ax_joint.tick_params(labelsize = fontsize)
ax_marg_x.tick_params(labelsize = fontsize)
ax_marg_y.tick_params(labelsize = fontsize)

pearsonr = str(np.round(stat.pearsonr(x, y)[0],decimals = 2))
pearsonp = str(np.int(stat.pearsonr(f_ksep, b_ksep)[1]))
ax_joint.text(-12, 26, 'pearsonr = ' + pearsonr + '; p = '+pearsonp, fontsize = fontsize)

plt.tight_layout()
plt.subplots_adjust(left = 0.135,
                    bottom = 0.11,
                    right = 0.99,
                    top = 0.99,
                    wspace = 0.11,
                    hspace = 0.11
)
plt.show()
fig.savefig(outdir + 'k1k2_Idif_modified.pdf', dpi = 1000)

