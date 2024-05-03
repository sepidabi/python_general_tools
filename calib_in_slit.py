from My_functions import *
import matplotlib.pyplot as ply
import scipy.io as sc
import math as math
import numpy as np
import lptools as lp
import satlas
import pyana as pyana
from scipy.interpolate import interp1d
from scipy.signal import fftconvolve
import crisp as fpi
from mpfit_my import mpfit
import sys



def fitme(par, x=True, y=True, imean=True, w=True, sig=True,fjac=None):
    fff= interp1d(x-par[1],y)
    y2= fff(w*par[2])*par[0]
    status = 0
    return [status,(y2-imean)*sig]

    
#choose the best frame
time_step=26


folder='/srv/scratch/crobu/sst/2017.04.20/'

line=[8542]

for i in line:

    if i == 8542: cw=8542.091
    else: cw=6302.4935
###################################################################################################
#load reference image

    fil0 = folder+'crispex.stokes.'+str(i)+'.09:40:25.time_corrected_CHROMIS.fcube'
    fil1 = folder+'crispex.stokes.'+str(i)+'.09:40:25.time_corrected_CHROMIS_sp.fcube'
    cube = read_polcube(fil0,fil1)
    nx=cube.shape[4]
    ny=cube.shape[3]
    nw=cube.shape[2]

#load f0 file and spectfile

    wl=(pyana.fzread(folder+'wav.'+str(i)+'.f0'))['data']
    spect_pos=sc.readsav(folder+'spectfile.'+str(i)+'.idlsave',verbose=True,python_dict=True)['spect_pos']

# load slit

    slit_file=folder+'slit_'+str(i)+'.csav'
    loop_slab=sc.readsav(slit_file,verbose=True,python_dict=True)['loop_slab']


###################################################################################################



#get solar atlas spectrum

    s = satlas.satlas()
    wav,inty,c = s.getatlas(wl[0]+cw-1.0,  wl[-1]+cw+1.0, cgs=True)




# select a QS patch from reference cube

    x0=911
    y0=68
    x1=1262
    y1=260

    qs_im=cube[time_step,0,:,y0:y1,x0:x1]



# Prepare the limb darkening correction

    mu= 0.33
    ptot=get_limbdarkening(mu, str(i))

    

#calculate the mean qs intensity for each wavelgt

    
    imean=np.zeros(nw)

    for ww in range (nw): imean[ww]=np.mean(qs_im[ww,:,:])

#convolve transmission profile with solar atlas        
        
    if i==8542:
        dum=np.argmin(inty)
        wav_delta= wav-wav[dum]
        ffpi=fpi.crisp(8542.0)
        dw=wav[1]-wav[0]
        numbp=51
        tw=(np.arange(numbp)-numbp/2)*dw 
        tr=ffpi.dual_fpi(tw, erh = -0.015)
        #tr /= tr.sum()
        inty_fpi = fftconvolve(inty*ptot, tr/np.sum(tr), mode='same')
        sig = np.arange(nw, dtype=np.float)+1.0    #to fix
        sig[4:12] = 0.3        #to fix

    if i==6302:
        int_index= ((wav>6302).nonzero())[0]
        dum=np.argmin(inty[int_index])
        wav_delta= wav-wav[int_index[dum]]
        ffpi=fpi.crisp(6302.0)
        dw=wav_delta[1]-wav_delta[0]
        numbp=51
        tw=(np.arange(numbp)-numbp/2)*dw 
        tr=ffpi.dual_fpi(tw, erh = -0.015)
        #tr /= tr.sum()
        inty_fpi = fftconvolve(inty*ptot, tr/np.sum(tr), mode='same')
        sig = np.ones(nw)
        sig[2:8] = 0.3 #to fix
        sig[11:16]=0.3
        sig[9] = 5   #to fix

    a={'fixed':0, 'limited':[1,1], 'limits':[0.,0.]}
    b={'fixed':0, 'limited':[1,1], 'limits':[0.,0.]}
    c={'fixed':0, 'limited':[1,1], 'limits':[0.,0.]}

    

    # a['limits']=[0., np.mean(imean)/np.mean(inty_fpi)*10.]
    # b['limits']= [-1.0,1.0]
    # c['limits']= [0.98,1.02]
    # c['fixed']= 1
    # fpar=[a,b,c]
    
    a['limits']=[0., np.mean(imean)/np.mean(inty_fpi[40:450])]
    b['limits']= [-0.5, 0.5]
    c['limits']= [1.0, 1.5]
    c['fixed'] = 1
    fpar=[a,b,c]

    
    
    farg = {'x':wav_delta[40:450], 'y':inty_fpi[40:450], 'w':wl, 'sig':sig, 'imean':imean}


    par= np.array([(np.mean(imean)/np.mean(inty_fpi[40:450])), 0.01, 1.0])
    print(par)

    #sys.exit()
    par = mpfit(fitme, par, parinfo = fpar, functkw = farg, debug=0)

    print 'status = ', par.status
    if (par.status <= 0): print 'error message = ', par.errmsg
    print 'parameters = ', par.params
    
    #par.params=[12653612.368581988,     0.015923775930778712,       1.0000000000000000]
    ply.plot(wav_delta, inty_fpi)
    ply.plot((wl+par.params[1])*par.params[2], imean/par.params[0])
    ply.show()

    


