import numpy as np
import lptools as lp
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from astropy.io import fits
from PIL import Image as image
import sparsetools as sp
import mfits as mf
import satlas
import crisp as fpi
import chromis as cr
from netCDF4 import Dataset as nc
from scipy.io import readsav as restore
import fnmatch
import os
import math as math
import imtools as im
import seaborn as sns
import scipy.stats as stat
import scipy.ndimage as spnd
import pickle
from ipdb import set_trace as stop
import scipy.ndimage as spnd
import scipy.signal as sgnl
from skimage import exposure
from scipy.fftpack import fft2, ifft2
import matplotlib.gridspec as gridspec


def get_name(var):
    for name, value in globals().items():
        if value is var:
            return name

def save_video():
    os.system("ffmpeg -r 1 -i img%01d.png -vcodec mpeg4 -y movie.mp4")

def file_search(filedir,pattern):
    files = os.listdir(filedir+'.')
    return sorted(fnmatch.filter(files,pattern))

#Carol + J (s-modification) recipe
def lp_read(fil0,fil1):

    nx, ny, dum, ns, dtype, ndim = lp.lphead(fil0)
    nw, nt, dum, ns, dtype, ndim = lp.lphead(fil1)

    io = np.memmap(fil0, mode='r', shape=(nt,ns,nw,ny,nx), dtype=dtype,
               offset=512)
    print '(nt, ns, nw, ny, nx) = (',nt,',',ns, ',',nw,',', ny,',',nx,')'
    return io


#def read_unpolcube(fil0,fil1):

#    nx, ny, dum, dtype, ndim = lp.lphead(fil0)
#    nw, nt, dum, dtype, ndim = lp.lphead(fil1)

#    io = np.memmap(fil0, mode='r', shape=(nt,nw,ny,nx), dtype=dtype,
#               offset=512)
#    print '(nt, nw, ny, nx) = (',nt,',',nw,',', ny,',',nx,')'
#    return io


def lp_read_scan(file):
    nx, ny, nt, dum, dtype, ndim = lp.lphead(file)
    io = np.memmap(file, mode='r', shape=(nt,ny,nx), dtype=dtype,
               offset=512)
    print '(nt, ny, nx) = (',nt,',', ny,',',nx,')'
    return io


def get_limbdarkening(mu, line):
    A={'6302':[0.33644, 1.30590, -1.79238, 2.45040, -1.89979, 0.59943, 0.8290],
       '6563':[0.34685, 1.37539, -2.04425, 2.70493, -1.94290, 0.55999, 0.8360],
       '8542':[0.45396, 1.25101, -2.02958, 2.75410, -2.02287, 0.59338, 0.8701],
       '3950':[0.12995, 0.91836, -0.07566, 0.19149, -0.28712, 0.12298, 0.7204]}
   
    ptot = 0.0

    for kk in range(6):
        ptot += (A[line][kk]*mu**(kk))

    return ptot


#def derivative(x,y):

#    n= x.shape[0]
 #   yp= np.zeros(n)
  #  
   # dx  = x[1:n-1] -x[0:n-2]
   # dx1 = x[2:] - x[1:n-1]
#
 #   der = (y[1:n-1] - y[0:n-2]) / dx
  #  der1 = (y[2:] - y[1:n-1]) / dx1

#    idx= np.where(der*der1 > 0.0)
 #   if (((der*der1 > 0.0).nonzero() >0) == True):
  #      lambdad = (1. + dx1[idx] / (dx1[idx] + dx[idx])) / 3.
   #     yp[np.asarray(idx)+1] = (der[idx] * der1[idx]) / (lambdad * der1[idx] + (1. - lambdad) * der[idx]) 
   # yp[0] = der[0]  
   # yp[n-1] = der1[der1.shape[0]-1]

#    return yp


def my_cluster(f0, f1, nclu, t_step=26, s=0, w0=0,w1=21, y0=50, x0=120):
    cube=read_polcube(f0,f1)
    nt=cube[:,s,w0:w1,y0:,x0:].shape[0]
    nw=cube[:,s,w0:w1,y0:,x0:].shape[1]
    ny=cube[:,s,w0:w1,y0:,x0:].shape[2]
    nx=cube[:,s,w0:w1,y0:,x0:].shape[3]
    mm= np.zeros((ny*nx))
    dat1=cube[t_step,s,w0:w1,y0:,x0:]
    dat =dat1.reshape(nw,nx*ny)
    mm=KMeans(n_clusters=nclu).fit(np.transpose(dat)).labels_
    return mm.reshape(ny,nx)


#Sepid's Recipe

#cgimage
def cgimage(im, ymin = 0, ymax = None, xmin = 0, xmax = None, xlabel = "[x]", ylabel = "[y]", cmap = "gray",title = "",show=False, save = False, savedir = '', name = 'image'): 

    if(xmax == None):
        xmax = im.shape[1]
    if(ymax ==None):
        ymax = im.shape[0]
        
    img = plt.imshow(im, cmap = cmap)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #img.set_tight_layout(True)
    plt.ion()
    plt.title(title)

    if(show == True):
        plt.show(img)

    if(save == True):
        plt.savefig(savedir+name)
        plt.close("all")
        print "saved to:  "+savedir+name

    return img

#cgplot
def cgplot(y, x = [],  ymin = None, ymax = None, xmin = None, xmax = None, xlabel = "[x]", ylabel = "[y]",title = "",show=False, color = 'black',sym = '', fill = True , symcolor = 'black', savedir = '', name = '',save = False,xtickstyle = '',ytickstyle = ''): 

 #   if(xmax == None):
  #      xmax = im.shape[1]
   # if(ymax ==None):
     #   ymax = im.shape[0]

    plot = plt.plot(x,y,sym, color = symcolor)
    if(fill == True):
        plot = plt.plot(x,y, color = color)
 
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ticklabel_format(style=xtickstyle,axis='x',)
    plt.ticklabel_format(style=xtickstyle,axis='y')
    #img.set_tight_layout(True)
    plt.ion()
    plt.title(title)
    if(show == True):
        plt.show()
    if(save==True):
         plt.savefig(savedir+name)
         plt.close('all')
         print "saved to:  "+savedir+name
    
    return plot

def readfits(nam,verbose = True):

    cub = fits.getdata(nam)
    print "dims = ", cub.shape
    if (verbose):
        print [fits.getheader(nam)]

    return cub

def unsharp(cube, sigma = 3., alpha = 1., mode = 'reflect', cmap = 'Greys', vmin = 0, vmax = 0, plot = False):
    cube_blur = spnd.gaussian_filter(cube, sigma, mode = mode)
    cube_unsharp = cube + alpha*(cube-cube_blur)
    if (vmin == 0 and vmax == 0):
        vmin, vmax = np.min(cube), np.max(cube)
    if(plot):
        plt.close('all')
        fig, (ax1, ax2) = plt.subplots(2, 1,figsize=(7,10))
        ax1.imshow(cube, cmap = cmap, vmin = vmin, vmax = vmax, origin = 'lower')
        ax2.imshow(cube_unsharp, cmap = cmap, vmin = vmin, vmax = vmax, origin = 'lower')
        plt.tight_layout()
        plt.show()
    return cube_unsharp

def derivative(xp, fxp, dx=1e-3):

    f2xp = np.zeros(len(xp))

    for i in range(len(xp)):
        if (i==0):
            f2xp[i] = (fxp[i+1] - fxp[i])/(xp[i+1] - xp[i])
        elif (i==len(xp)-1):
            f2xp[i] = (fxp[i] - fxp[i-1])/(xp[i] - xp[i-1])
        else:
            f2xp[i] = (fxp[i+1] - fxp[i-1])/(xp[i+1] - xp[i-1])

    x = np.linspace(xp[0], xp[-1], int((xp[-1]-xp[0])/dx))
    f2x = np.interp(x,xp,f2xp)
    fx = np.interp(x,xp,fxp)

    return x, fx, f2x

def integral(xp, fxp, dx=1e-3):

    f2xp = np.zeros(len(xp))

    for i in range(len(xp)):
        if (i==0):
            f2xp[i] = (fxp[i+1] + fxp[i])*(xp[i+1] - xp[i])/2.
        elif (i==1):
            f2xp[i] = f2xp[i-1] + fxp[i]*(xp[i] - xp[i-1])
        else:
            f2xp[i] = f2xp[i-2] + fxp[i-1]*(xp[i] - xp[i-2])

    x = np.linspace(xp[0], xp[-1], int((xp[-1]-xp[0])/dx))
    f2x = np.interp(x,xp,f2xp)
    fx = np.interp(x,xp,fxp)

    return x, fx, f2x


def gamma(im,g,c=10):
    im_corr = np.e**(g*np.log(im)+np.log(c))
    return im_corr

def save_obj(filename, obj, mode = "wb"):
    pickling_on = open(filename, mode)
    pickle.dump(obj, pickling_on)
    pickling_on.close()
    print("object saved to: " + filename)

def load_obj(filename, mode = "rb"):
    pickle_off = open(filename, mode)
    obj = pickle.load(pickle_off)

    return obj

def image_int_correct(cube, mode, kernel = [1,1], sigma=None, contrast ='log', gamma_f = 1, log_f = 1, show = False):
    
    if mode=='median':
        cube_correct = sgnl.medfilt2d(cube,kernel)

    if mode=='gauss':
        if sigma==None:
            sigma = np.std(cube)
        cube_correct = spnd.gaussian_filter(cube,sigma)

    if mode=='hist_eq':
        cube_correct = exposure.equalize_hist(cube)

    if contrast=='gamma':
        cube_final = exposure.adjust_gamma(cube_correct,gamma_f)

    if contrast=='log':
        cube_final = exposure.adjust_log(cube_correct,log_f)

    if (show):
        plt.close('all')
        test = plt.figure(figsize = (15,8))
        gs = gridspec.GridSpec(1,2, figure = test)
        ax_before = plt.subplot(gs[0,0])
        ax_after = plt.subplot(gs[0,1])

        ax_before.imshow(cube, cmap = 'gray', origin = 'lower')
        ax_after.imshow(cube_final,cmap = 'gray', origin = 'lower')
        plt.tight_layout()
        plt.show()

    return cube_final

def sharpen(f, sigma=1, unsigma=1, alpha = 30):
    blurred_f = spnd.gaussian_filter(f, sigma = sigma)

    filter_blurred_f = spnd.gaussian_filter(blurred_f, sigma = unsigma)

    alpha = 30
    sharpened = blurred_f + alpha * (blurred_f - filter_blurred_f)
    return sharpened

