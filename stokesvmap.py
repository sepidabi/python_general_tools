from My_functions import *
import numpy as np
import scipy.io as sc
import matplotlib
matplotlib.use("TkAgg")
import Tkinter as Tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import math as math
import scipy as sci
import idlsave
import matplotlib.pyplot as plt
import sparsetools as sp
import imtools as im
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
from ipdb import set_trace as stop
from matplotlib.widgets import MultiCursor



class iGui:
    def __init__(self, mapp, v,nw):
        slider_length = 250
        self.mapp = mapp
        self.v = v

        self.x = 200
        self.y = 100
        

        plt.ion()

        self.fig1, self.ax1 = plt.subplots(figsize=(8,5), nrows = 1, ncols = 1, dpi = 100)
        self.fig2, self.ax2 = plt.subplots(figsize=(8,5), nrows = 1, ncols = 1, dpi = 100)


        self.im01 = [None]
        self.im01  = self.ax1.imshow(self.mapp, cmap = 'seismic', origin = 'lower', vmin=-20,vmax=20)
        self.im02, = [None]
        self.im02, = self.ax2.plot(np.arange(nw),v[:,self.y,self.x].squeeze())
        self.ax2.set_ylim(-100,100)

        
        self.cid = self.fig1.canvas.mpl_connect('button_press_event', self.getXY)
        self.multi = MultiCursor(self.fig1.canvas, (self.ax1, self.ax2), color='w', lw=0.75,useblit=True, vertOn=True, horizOn=True)

        self.fig1.show()
        self.fig1.set_tight_layout(True)
        self.fig2.show()
        self.fig2.set_tight_layout(True)


    def getXY(self, event):
        if(event.inaxes == self.ax1):
            self.x = int(event.xdata)
            self.y = int(event.ydata)
            print("(x,y) = ({0},{1})".format(self.x, self.y))
            self.reDrawPlots()

        
        
    def reDrawPlots(self):
        print(self.v[:,self.y,self.x].squeeze())
        self.im02.set_data(np.arange(nw),self.v[:,self.y,self.x])
        self.fig2.canvas.draw()
        self.fig2.canvas.flush_events()

        
if __name__ == "__main__":


    dir='/srv/scratch/crobu/sst/2017.04.20/'
    file0=dir+'crispex.stokes.8542.09:40:25.time_corrected_CHROMIS.fcube'
    file1=dir+'crispex.stokes.8542.09:40:25.time_corrected_CHROMIS_sp.fcube'
    spect_pos=sc.readsav(dir+'spectfile.8542.idlsave',verbose=True,python_dict=True)['spect_pos']
    
    # file0=dir+'crispex.stokes.6302.09:40:25.time_corrected_CHROMIS.fcube'
    # file1=dir+'crispex.stokes.6302.09:40:25.time_corrected_CHROMIS_sp.fcube'
    # spect_pos=sc.readsav(dir+'spectfile.6302.idlsave',verbose=True,python_dict=True)['spect_pos']

    nw=spect_pos.shape[0]
    
    # tsteps=np.append(np.append(np.arange(5,14),np.arange(24,28)),np.append(np.arange(32,38),np.arange(42,46)))

    cube=read_polcube(file0,file1)
    
    #V=np.median(cube[tsteps,3,:,:,:],axis=0)
    # V_b = np.sum(cube[:,3,0:6,:,:], axis=0)
    # V_r = np.sum(cube[:,3,6:,:,:], axis=0)

    v_tot=np.sum(cube[:,3,:,:,:], axis=0)
    
    
    # mapp = np.zeros((v_tot.shape[1],v_tot.shape[2])

    # for y in range (V_b.shape[1]):
    #     print(y,V_b.shape[1]-1)
    #     for x in range (V_b.shape[2]):    
    #                  if (v_tot[6,y,x] > 0): mapp[y,x]=1

        
    #v_tot=np.append(V_b,V_r, axis=0)
    
    
    #app = iGui(cube[26,3,6,:,:],v_tot,nw)
    plt.imshow(v_tot[6,:,300:], cmap='gist_gray', origin='lower')
    plt.contour(v_tot[6,:,300:],)
    plt.show()
    #app = iGui(mapp,v_tot,nw)
    

 
