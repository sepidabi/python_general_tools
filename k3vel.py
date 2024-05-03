from examply import read_polcube
import matplotlib.pyplot as ply
import scipy.io as sc
from scipy.optimize import basinhopping
import math as math
import numpy as np
import lptools as lp
from scipy.constants import *
from scipy import interpolate
from astropy.modeling import models, fitting


#files
folder='/srv/scratch/crobu/sst/2017.04.20/'
fil0 = folder+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected.fcube'
fil1 = folder+'crispex_3950_2017-04-20T09:41:33_scans=0-75_time-corrected_sp.fcube'


cube = read_polcube(fil0,fil1)

print("cube (nt,nw,ny,nx)=",cube.shape)

nx=cube.shape[4]
ny=cube.shape[3]
nw=cube.shape[2]
nt=cube.shape[0]

spect_pos=sc.readsav(folder+"spectfile.3950.idlsave",verbose=True,python_dict=True)['spect_pos']


minimizer_kwargs = {"method": "BFGS"}

tt= 4
wrange=np.arange(4,15)
spect_pos=spect_pos[wrange]
x0=[7.]
s=(ny,nx)
minimo=np.zeros(s)
vel=np.zeros(s)
interp_array=np.arange(1000)/999.*(spect_pos[-1]-spect_pos[0])+spect_pos[0]


for xx in range(10):
    print(10)
    for yy in range(ny):
        
        r=basinhopping(func, x0)
        
        g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, spect_pos, cube[tt,0,wrange,xx,yy])
        interpf=interpolate.interp1d(spect_pos,g)
        new_g=interpf(interp_array)
        minimo[xx,yy]= (np.abs(new_g[:] - min(new_g[:])) < 1e-14).nonzero()

        

#         vel[yy,xx]=(1- spect_pos[9]/spect_pos[minimo[yy,xx]])*c/1000. #km

# ply.plot(cube[tt,0,0:20,yy,xx])
# ply.plot([minimo, minimo], [1e-9,cube[tt,0,20,yy,xx]], color='k', linestyle='-', linewidth=2)

# ply.show()
# print(minimo)

# ply.imshow(vel)
# ply.show()
