#CHI2_TEST program designed to check the fitting of
#different inversion settings

import sparsetools as sp
import matplotlib.pyplot as plt
import numpy as np
from sepid import *
#import spectral as s

#DECLERATIONS
index = '06_165307'
index_bg = '06_165404'
ss = 0 # stokes param index
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
outdir = '/home/seki2695/OUTPUT/inv/'


reg = np.array([0.1,0.316,1.,3.16,10.0],dtype = float)
reg_type = ['1,2,3,1,1,0,1']*len(reg)
reg_w = ['1,1,1,1,1,1,1']*len(reg)
chi2_fm = np.zeros(len(reg))
chi2_fx = np.zeros(len(reg))
chi2_fmed = np.zeros(len(reg))
chi2_fs = np.zeros(len(reg))

mod = ['VIII', 'V', 'VI', 'VII', 'IV']
cyc = ['1','1','1','1','1']
amp = 1. #cak amplifier

for m in range(len(reg)):
    plt.close('all')
    reg_ann = 'n$_{T}$ = -6.8, -5.9, -5.2, -4.3, -3.4, -2.4, -1.5, -0.8, 0.08 \nn$_{v_{LOS}}$ =  -6.5, -5.0, -4.0, -3.0, -1.7, -0.2 \nn$_{v_{turb}}$ = -6.0, -4, -1.9, -0.2 \nreg$_{type}$ = '+reg_type[m]+'\nreg =' +"%.2f" % round(reg[m],2)+'\nreg$_{weight}$ = '+reg_w[m]+'\n$\mathrm{\chi^2_{f}}$ = '

    #Extracting inversion results
    f_i = sp.profile('fb/f'+mod[m]+'_observed_'+index+'.nc')
    f_o = sp.profile('fb/f'+mod[m]+'_synthetic_cycle'+cyc[m]+'_'+index+'.nc')
    #f_m = sp.model('fb/f'+mod[m]+'_atmosout_cycle'+cyc[m]+'_'+index+'.nc')
    b_i = sp.profile('bg/b'+mod[m]+'_observed_'+index+'.nc')
    b_o = sp.profile('bg/b'+mod[m]+'_synthetic_cycle'+cyc[m]+'_'+index+'.nc')
    #b_m = sp.model('bg/b'+mod[m]+'_atmosout_cycle'+cyc[m]+'_'+index+'.nc')

    #Fibril extract
    #========
    #inversion
    f_obs = f_i.dat[0,0,:,:,ss]
    f_obs[:,0:42]  = f_i.dat[0,0,:,0:42,ss]*amp
    f_syn = f_o.dat[0,0,:,:,ss]
    f_syn[:,0:42]  = f_o.dat[0,0,:,0:42,ss]*amp
    f_wav = f_i.wav
    f_n = f_i.dat.shape[2]
    #f_dep = f_m.ltau[0,0,0,:]
    #slab
    #fibdir = datadir+'fr'+str(fn)+'/'
    #fib_file = file_search(fibdir,'crispex*3950*'+index+'*.csav')
    #fib = restore(fibdir+fib_file[0])
    #f_slab = fib.loop_slab  #intensity values in the desired frame
    #f_x = fib.x_coords
    #f_x_pts = fib.x_loop_pts
    #f_y = fib.y_coords
    #f_l = fib.loop_size
    #f_y_pts = fib.y_loop_pts

    #Background extract
    #===========
    #inversion
    b_obs = b_i.dat[0,0,:,:,ss]
    b_obs[:,0:42] = b_i.dat[0,0,:,0:42,ss]*amp
    b_syn = b_o.dat[0,0,:,:,ss]
    b_syn[:,0:42] = b_o.dat[0,0,:,0:42,ss]*amp
    b_wav = b_i.wav
    b_n = b_i.dat.shape[2]
    #b_dep = b_m.ltau[0,0,0,:]
    #slab
    #bg_file = file_search(fibdir,'crispex*3950*'+index_bg+'*.csav')
    #bg = restore(fibdir+bg_file[0])
    #b_slab = bg.loop_slab #intensity values in the desired frame
    #b_x = bg.x_coords
    #b_x_pts = bg.x_loop_pts
    #b_y_pts = bg.y_loop_pts
    #b_y = bg.y_coords
    #b_l = bg.loop_size

    #chi2
    pxnl = np.min([f_n,b_n])
    chi2_f = np.zeros(pxnl)
    chi2_b = np.zeros(pxnl)
    for px in range(0,pxnl):
        chi2_f[px] = np.sum(((f_o.dat[0,0,px,:,ss] - f_i.dat[0,0,px,:,ss])/f_o.weights[:,ss])**2)/len(f_wav)
        chi2_b[px] = np.sum(((b_o.dat[0,0,px,:,ss] - b_i.dat[0,0,px,:,ss])/b_o.weights[:,ss])**2)/len(f_wav)
    chi2_fm[m] = np.mean(chi2_f)
    chi2_fx[m] = np.amax(chi2_f)
    chi2_fmed[m] = np.median(chi2_f)
    chi2_fs[m] = np.std(chi2_f)
    
    plt.plot(chi2_f, color = 'red', label = 'Fibril')
    plt.plot(chi2_b, color = 'blue', label = 'Background')
    plt.xlabel(r'N$_\mathrm{data}$')
    plt.ylabel(r'$\chi^2$')
    plt.ylim([0,30])
    leg = plt.legend()
    plt.title(reg_ann + "%.2f" % round(np.mean(chi2_f),2))
    plt.tight_layout()
    #bb = leg.legendPatch.get_bbox().inverse_transformed(plt.transAxes)
    #plt.text(bb.x0, bb.y0 - bb.height, 'text', transform=plt.transAxes)
    #plt.tex()

    plt.savefig(outdir+'chi2_test_mod'+mod[m]+'.pdf')
    print 'file saved to: '+ outdir+'chi2_test_mod'+mod[m]+'.pdf'
    #plt.show()

plt.close('all')
plt.plot(np.log10(np.sqrt(reg)),chi2_fm, label = r'$\mathrm{\overline{\chi^2}}$')
plt.plot(np.log10(np.sqrt(reg)),chi2_fm,'.')
plt.plot(np.log10(np.sqrt(reg)),chi2_fx, label = r'$\mathrm{max\{\chi^2\}}$')
plt.plot(np.log10(np.sqrt(reg)),chi2_fx,'.')
#plt.plot(np.log10(np.sqrt(reg)),chi2_fmed, label = r'$\mathrm{[\chi^2]_{med}}$')
#plt.plot(np.log10(np.sqrt(reg)),chi2_fmed,'.')
#plt.plot(np.log10(np.sqrt(reg)),chi2_fs, label = r'$\mathrm{\sigma_{\chi^2}}$')
#plt.plot(np.log10(np.sqrt(reg)),chi2_fs,'.')
plt.legend()
plt.xlabel(r'$\mathrm{log(\sqrt{\lambda^2})}$')
plt.xticks([-0.5,-0.25,0,0.25,0.5])
plt.xlim([-0.55,0.55])
plt.ylabel(r'$\mathrm{\chi^2}$')
plt.savefig(outdir+'chi2_stat.pdf')
