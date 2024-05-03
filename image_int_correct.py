# fibril_hor_oscillation.py
'''meant to improve the image intensity
to bring out the features,
e.g. the bright fibrils better
'''

# import modules
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
import numpy as np
import scipy as sci
import scipy.ndimage as spnd
import scipy.signal as sgnl
from skimage import exposure
from scipy.fftpack import fft2, ifft2
import matplotlib.gridspec as gridspec

def image_int_correct(cube, mode, kernel = [1,1], gauss_sigma=None, contrast ='log', gamma_f = 1, log_f = 1, show = False):
    
    if mode=='median':
        cube_correct = sgnl.medfilt2d(cube,kernel)

    if mode=='gauss':
        if gauss_sigma==None:
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

        

        
        
