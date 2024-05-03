from My_functions import *
import numpy as np
import matplotlib.pyplot as ply
import scipy.io as sc
import math as math
import lptools as lp
import satlas
import pyana as pyana
from scipy.interpolate import interp1d
from scipy.signal import fftconvolve
import crisp as fpi
from mpfit_my import mpfit
import sys
import imtools as im
import os
from astropy.io import fits
import sparsetools as sp
from netCDF4 import Dataset as nc
import chromis as cr
from ipdb import set_trace as stop
import mathtools as mt

#-----------------------------------------------------------
def mfits(file):
    a = fits.open(file, 'readonly')
    return np.float32(a[0].data[:])

#-----------------------------------------------------------

def findgrid(w, dw, odw,extra):
    w1=np.round(w*1000).astype('int32')
    dw1 = int(dw*1000)

    w2 = w/dw
    w2 = np.round(w2)
    
    idx = np.zeros(w.size, dtype='int32')

    np0 = w2[-1] - w2[0] + extra +1
    wn = np.arange(np0, dtype='float64')*odw + w[0] -extra/2*odw

    for ii in range(w.size):
        idx[ii] = np.argmin(np.abs(wn-w[ii]))

    

    # ingrid = np.int32(np.round(w / dw))
    # w = ingrid * dw

    # ma = w.max()
    # mi = w.min()
    # n = int(np.round((ma-mi)/odw)) + 1 + extra
    
    # wn = np.arange(n, dtype='float64')*odw + mi - extra/2*odw
    # idx = np.zeros(w.size, dtype='int32')
    # for ii in range(w.size):
    #     idx[ii] = np.argmin(np.abs(wn-w[ii]))
    
    return wn, idx

#-----------------------------------------------------------    

def writeInstProf(oname, var, pref=[]):
    ncfile1 = nc(oname,'w', format='NETCDF4')
    ncfile1.createDimension('wav',var.size)
    par1 = ncfile1.createVariable('iprof','f8',('wav'))
    par1[:] = var


    if(len(pref) == 3):
        ncfile1.createDimension('np',len(pref))
        par2 = ncfile1.createVariable('pref','f8',('np'))
        par2[:] = np.float32(pref)

    ncfile1.close()

#-----------------------------------------------------------
def getCont(lam):
    s = satlas.satlas()
    x,y,c = s.getatlas(lam-0.1,lam+0.1, cgs=True)
    return np.median(c)
#-----------------------------------------------------------


#choose the best frame

cw_ca8=8542.091
cw_fe=6302.4935
cw_cak=3933.6640
cw_cah=3968.4690
time_step=26
mu=0.33
cont = {'cak':getCont(cw_cak)/get_limbdarkening(mu,'3950'), 'cah':getCont(cw_cah)/get_limbdarkening(mu,'3950'),
        'ca8':getCont(cw_ca8)/get_limbdarkening(mu,'8542'), 'fe':getCont(cw_fe)/get_limbdarkening(mu,'6302')}

folder='/srv/scratch/crobu/sst/2017.04.20/'
odir='/srv/scratch/crobu/sst/2017.04.20/prep_4stic/'
 

#load reference image

fil0_ca8 = folder+'crispex.stokes.8542.09:40:25.time_corrected_CHROMIS.fcube'
fil1_ca8= folder+'crispex.stokes.8542.09:40:25.time_corrected_CHROMIS_sp.fcube'
cube_ca8 = read_polcube(fil0_ca8,fil1_ca8)
nx=cube_ca8.shape[4]
ny=cube_ca8.shape[3]
nw=cube_ca8.shape[2]

fil0_fe = folder+'crispex.stokes.6302.09:40:25.time_corrected_CHROMIS.fcube'
fil1_fe = folder+'crispex.stokes.6302.09:40:25.time_corrected_CHROMIS_sp.fcube'
cube_fe = read_polcube(fil0_fe,fil1_fe)

fil0_cahk = folder+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected.fcube'
fil1_cahk = folder+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected_sp.fcube'
cube_cahk = read_polcube(fil0_cahk,fil1_cahk)

#load f0 file and spectfile

wl_ca8=(pyana.fzread(folder+'wav.8542.f0'))['data']
spect_pos_ca8=sc.readsav(folder+'spectfile.8542.idlsave',verbose=True,python_dict=True)['spect_pos']

wl_fe=(pyana.fzread(folder+'wav.6302.f0'))['data']
spect_pos_fe=sc.readsav(folder+'spectfile.6302.idlsave',verbose=True,python_dict=True)['spect_pos']

wl_cahk=((pyana.fzread(folder+'wav.3950.f0'))['data'])*1e10 #in angstrom
spect_pos_cahk=(sc.readsav(folder+'spectfile.3950.idlsave',verbose=True,python_dict=True)['spect_pos'])

# load slit

slit_file_ca8=folder+'slit_8542.csav'
loop_slab_ca8=sc.readsav(slit_file_ca8,verbose=True,python_dict=True)['loop_slab']

slit_file_fe=folder+'slit_6302.csav'
loop_slab_fe=sc.readsav(slit_file_fe,verbose=True,python_dict=True)['loop_slab']
xp = sc.readsav(slit_file_fe,verbose=True,python_dict=True)['x_loop_pts']
yp = sc.readsav(slit_file_fe,verbose=True,python_dict=True)['y_loop_pts']


#find grid

wfe1, ife1 = findgrid(wl_fe[0:10], 0.01, 0.01, extra=8)
wfe2, ife2 = findgrid(wl_fe[10:], 0.01, 0.01, extra=12)
wfe = np.append(wfe1,wfe2)
ife   = np.append(ife1, ife2+wfe1.size)
wca8, ica8= findgrid(wl_ca8, 0.05,0.05, extra=8)
wck, ick = findgrid(wl_cahk[0:21], (wl_cahk[10]-wl_cahk[9]),(wl_cahk[10]-wl_cahk[9]), extra=8)
wch, ich = findgrid(wl_cahk[21:-1], (wl_cahk[31]-wl_cahk[30]),(wl_cahk[31]-wl_cahk[30]), extra=8)



#load wav_calibration files

calib_8542 = mfits('calib.8542.fits')
calib_6302 = mfits('calib.6302.fits')
calib_3934 =  mfits('calib.3934.fits')
calib_3968 =  mfits('calib.3968.fits')

wlcalib_8542  =(wca8+calib_8542[1])*calib_8542[2]+cw_ca8
wlcalib_6302_1=(wfe1+calib_6302[1])*calib_6302[2]+cw_fe
wlcalib_6302_2=(wfe2+calib_6302[1])*calib_6302[2]+cw_fe
wlcalib_3934 = (wck+ calib_3934[1])*calib_3934[2]
wlcalib_3968 = (wch+ calib_3968[1])*calib_3968[2]
wlcalib_cont = (wl_cahk[-1]+ calib_3968[1])*calib_3968[2]

    
fe_1   = sp.profile(nx = loop_slab_fe.shape[2],  ny=1, ns=4, nw=wfe1.size)
fe_2   = sp.profile(nx = loop_slab_fe.shape[2],  ny=1, ns=4, nw=wfe2.size)
ca_8   = sp.profile(nx = loop_slab_ca8.shape[2], ny=1, ns=4, nw=wca8.size)
ca_k   = sp.profile(nx = loop_slab_ca8.shape[2], ny=1, ns=4, nw=wck.size)
ca_h   = sp.profile(nx = loop_slab_ca8.shape[2], ny=1, ns=4, nw=wch.size+1)


fe_1.wav[:] = wlcalib_6302_1[:]
fe_2.wav[:] = wlcalib_6302_2[:]#+0.02
ca_8.wav[:] = wlcalib_8542[:]-0.06
ca_k.wav[:] =  wlcalib_3934[:]  #wck[:]
ca_h.wav[0:-1]= wlcalib_3968[:] #wch[:]
ca_h.wav[-1]  = wlcalib_cont   #wl_cahk[-1] # Continuum point

for ii in range(4):
    fe_1.dat[0,0,:,ife1,ii] =  np.transpose(cube_fe[time_step,ii,0:10, yp.astype(int), xp.astype(int)] / calib_6302[0]/cont['fe'])
    fe_2.dat[0,0,:,ife2,ii] =  np.transpose(cube_fe[time_step,ii,10:, yp.astype(int), xp.astype(int)] / calib_6302[0]/cont['fe'])
    ca_8.dat[0,0,:,ica8,ii] = np.transpose(cube_ca8[time_step,ii,:, yp.astype(int), xp.astype(int)]/ calib_8542[0]/cont['ca8'])

ca_k.dat[0,0,:,ick[:],0] = np.transpose(cube_cahk[time_step,0,0:21, yp.astype(int), xp.astype(int)])/calib_3934[0] / cont['cak']
ca_h.dat[0,0,:,ich[:],0] = np.transpose(cube_cahk[time_step,0,21:-1, yp.astype(int), xp.astype(int)])/calib_3968[0] / cont['cah'] 
ca_h.dat[0,0,:,-1,0] = np.transpose(cube_cahk[time_step,0,-1, yp.astype(int), xp.astype(int)]) / getCont(wl_cahk[-1])


fe_1.weights[:,:] = 1.e32 # Very high value means weight zero
fe_1.weights[ife1,:] = 0.005
fe_1.weights[ife1,1:3] /= 1. # Some more weight for Q&U
fe_1.weights[ife1,3] /= 1.    # Some more weight for V

fe_2.weights[:,:] = 1.e32 # Very high value means weight zero
fe_2.weights[ife2,:] = 0.005
fe_2.weights[ife2,1:3] /= 1. # Some more weight for Q&U
fe_2.weights[ife2,3] /= 1.   # Some more weight for V
fe_2.weights[ife2[-3::],:] = 1.e32 # Telluric



ca_8.weights[:,:] = 1.e32 # Very high value means weight zero
ca_8.weights[ica8,:] = 0.005
ca_8.weights[ica8,1:3] /= 1.0 # Some more weight for Q&U
ca_8.weights[ica8,3] /= 1.0    # Some more weight for V


ca_k.weights[:,:] = 1.e32 # Very high value means weight zero
ca_k.weights[ick,0] = 0.003
#ca_k.weights[ick[7:12],0] /= 2.0

ca_h.weights[:,:] = 1.e32 # Very high value means weight zero
ca_h.weights[ich,0] = 0.003
#ca_h.weights[ich[7:12],0] /= 2.0
ca_h.weights[-1,0] = 0.004 # Continuum point


#sp_all = ca_k + ca_h + fe_1 + fe_2 + ca_8
sp_all = fe_1 + fe_2 + ca_8
sp_all.write(odir+'observed_whk.nc')
#sp_all.write(odir+'observed.nc')
#sp_altern= sp_all.skip(6,1)
#sp_altern.write(odir+'observed_reduced.nc')

kk = sp_all.extractPix(x0=200, x1=201, y0=0, y1=1)
kk.write(odir+'observed_1pix.nc')

kk.dat[0,0,0,:,:] = sp_all.averageSpectrum()
kk.write(odir+'observed_average.nc')


# Ca II 8542
dw =  ca_8.wav[1]-ca_8.wav[0]
ntw= 25
f=fpi.crisp(8542.0)
tw = (np.arange(ntw)-ntw/2)*dw 
tr = f.dual_fpi(tw, erh = -0.01)
tr /= tr.sum()
writeInstProf(odir+'8542.nc', tr,  [8542.091, 9.0, 2.0])

    
# 6301/6302, we will use the same profile, so it should not have more points than any
# of the 2 regions
dw =  fe_1.wav[1]-fe_1.wav[0]
ntw= 45 
f=fpi.crisp(6302.0)
tw = (np.arange(ntw)-ntw/2)*dw 
tr = f.dual_fpi(tw, erh=-0.01)
tr /= tr.sum()
writeInstProf(odir+'6302.nc', tr,  [6302.1, 4.4, 2.0])


# Ca K region
dw =  ca_k.wav[3]-ca_k.wav[2]
ntw= 21 # Always an odd number < ca_k.wav.size
tw1 = (np.arange(ntw)-ntw/2)*dw + cw_cak
tr1 = cr.dual_fpi(tw1, erh = -0.1)
tr1 /= tr1.sum()
# Stores the FPI profile and the parameters of the prefilter
writeInstProf('3934.nc', tr1, [ca_k.wav[10], 4.5, 3.0])


# Ca H region
dw =  ca_h.wav[3]-ca_h.wav[2]
ntw= 21 # Always an odd number < ca_k.wav.size
tw1 = (np.arange(ntw)-ntw/2)*dw + cw_cah
tr1 = cr.dual_fpi(tw1, erh = -0.1)
tr1 /= tr1.sum()
# Stores the FPI profile and the parameters of the prefilter
writeInstProf('3968.nc', tr1, [ca_h.wav[10], 4.5, 3.0])

    

# First create a tau scale

mi = -7.6
ma = 1.0
dt = 0.18
ntau = int((ma-mi)/dt+0.5) + 1
tau = np.arange(ntau)/(ntau-1.0) * (ma-mi) + mi
itau = np.float64((-7.6, -5.0, -3.0, -1.0, 1.0))
itemp = np.float64((25., 5.5, 4.0, 4.5, 7.0 ))*1000
temp = mt.bezier3(itau, itemp, tau)   


# Fill in the model
m = sp.model(nx=loop_slab_fe.shape[2],  ny=1, nt=1, ndep=ntau)



# The inversion only needs to know the gas pressure at the upper boundary. FALC has Pgas[top] ~ 0.3, but
# this value is for quiet-Sun. Active regions can have up to Pgas[top] = 10.

m.pgas[:,:,:,0] = 1.0    

# Fill in initial B field and velovity (optional)
m.vturb[:,:,:,:] = 1.e5
m.vlos[:,:,:,:] = 0.5e5 # cm/s
m.Bln[:,:,:,:] = 400.
m.Bho[:,:,:,:] = 400.
m.azi[:,:,:,:] = 45. * 3.14159 / 180.

for ii in range(ntau):
    m.ltau[:,:,:,ii] = tau[ii]
    m.temp[:,:,:,ii] = temp[ii]

# Write to HD
m.write(odir+'modelin.nc') #, write_all=True)
m1 = m.extract(x0=0, x1=1, y0=0, y1=1)
m1.write(odir+'modelin_1pix.nc')

#m1 = m.scale(nx=sp_altern.nx, ny=sp_altern.ny)
#m1.write(odir+'modelin_reduced.nc')


#
# Now print the regions for the config file
#
lab = "region = {0:10.5f}, {1:8.5f}, {2:3d}, {3:e}, {4}"
print(" ")
print("Regions information for the input file:" )
#print(lab.format(ca_k.wav[0], ca_k.wav[1]-ca_k.wav[0], ca_k.wav.size, cont['cak'], 'fpi, 3934.nc'))
#print(lab.format(ca_h.wav[0], ca_h.wav[1]-ca_h.wav[0], ca_h.wav.size-1, cont['cah'], 'fpi, 3968.nc'))
#print(lab.format(ca_h.wav[-1], ca_h.wav[1]-ca_h.wav[0], 1, getCont(wl_cahk[-1]), 'none, none'))
print(lab.format(fe_1.wav[0], fe_1.wav[1]-fe_1.wav[0], fe_1.wav.size, cont['fe'], 'fpi, 6302.nc'))
print(lab.format(fe_2.wav[0], fe_2.wav[1]-fe_2.wav[0], fe_2.wav.size, cont['fe'], 'fpi, 6302.nc'))
print(lab.format(ca_8.wav[0], ca_8.wav[1]-ca_8.wav[0], ca_8.wav.size, cont['ca8'], 'fpi, 8542.nc'))
print("(w0, dw, nw, normalization, degradation_type, instrumental_profile file)")
print(" ")
    

f = ply.figure()
ax1 = ply.subplot2grid((2,2), (0,0))
ax2 = ply.subplot2grid((2,2), (0,1))
ax3 = ply.subplot2grid((2,2), (1,0))
ax4 = ply.subplot2grid((2,2), (1,1))
#ax1.set_yticklabels([])
#ax1.set_xticklabels([])
#ax2.set_yticklabels([])
#ax2.set_xticklabels([])
#ax3.set_yticklabels([])
#ax3.set_xticklabels([])
#ax4.set_yticklabels([])
#ax4.set_xticklabels([])

pt=200

ax1.plot(fe_1.wav[ife1], fe_1.dat[0,0,pt, ife1,0],'ro')
fnct=interp1d(fe_1.wav[ife1], fe_1.dat[0,0,pt, ife1,0],fill_value='extrapolate')
ax1.plot(fe_1.wav[:], fnct(fe_1.wav[:]),'b')

ax1.plot(fe_2.wav[ife2], fe_2.dat[0,0,pt, ife2,0],'ro')
fnct=interp1d(fe_2.wav[ife2], fe_2.dat[0,0,pt, ife2,0],fill_value='extrapolate')
ax1.plot(fe_2.wav[:], fnct(fe_2.wav[:]),'b')

ax2.plot(ca_8.wav[ica8], ca_8.dat[0,0,pt, ica8,0],'ro') 
fnct=interp1d(ca_8.wav[ica8], ca_8.dat[0,0,pt, ica8,0],fill_value='extrapolate')
ax2.plot(ca_8.wav[:], fnct(ca_8.wav[:]),'b')

ax3.plot(ca_k.wav[ick[0:-1]], ca_k.dat[0,0,pt, ick[0:-1],0],'ro')
fnct=interp1d(ca_k.wav[ick[0:-1]], ca_k.dat[0,0,pt, ick[0:-1],0],fill_value='extrapolate')
ax3.plot(ca_k.wav[0:-1], fnct(ca_k.wav[0:-1]),'b')

ax4.plot(ca_h.wav[ich[0:-1]], ca_h.dat[0,0,pt, ich[0:-1],0],'ro')
fnct=interp1d(ca_h.wav[ich[0:-1]], ca_h.dat[0,0,pt, ich[0:-1],0],fill_value='extrapolate')
ax4.plot(ca_h.wav[0:ich[-2]], fnct(ca_h.wav[0:ich[-2]]),'b')



ply.show()
