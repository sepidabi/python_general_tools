#!/usr/bin/env python 
import matplotlib
matplotlib.use("TkAgg")

import Tkinter as Tk
import math
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import sparsetools as sp
import imtools as im
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
from ipdb import set_trace as stop
from matplotlib.widgets import MultiCursor


plt.rcParams['toolbar'] = 'None'

class iGui:
    def __init__(self, root, obs, syn, mod, t0 = 0):
        slider_length = 250
        # read data
        self.m = sp.model(mod)
        self.s = sp.profile(syn)
        self.o = sp.profile(obs)
#        self.m. *= 180./3.1415926
        self.m.azi *= 180./3.1415926
        self.m.vturb *= 1.e-5
        self.m.vlos *= 1.e-5
        self.t0 = t0
        
        # Init
        self.xx = 0
        self.yy = 0
        self.idx = np.where(self.o.dat[0,0,0,:,0] > 0.0)
        self.w = np.arange(self.o.nw)
        self.itau = np.argmin(np.abs(self.m.ltau[0,0,0]))
        m = self.m
        self.vars = [m.temp, m.vlos, m.vturb, m.Bln, m.Bho, m.azi, np.log10(m.pgas), np.log(m.nne)]
        self.names= ['Temp','Vlos','Micro', 'LnB', 'HoB', 'Azi','Pgas','Ne']
        self.v1 = self.vars[0]
        self.v2 = self.vars[1]
        self.thres = 1.e-3
        
        # Init menu
        self.master = root
        self.master.wm_title("Inspector")
        self.frame = Tk.Frame(self.master)  
        
        # create Figures
        plt.ion()
        self.fig, self.ax = plt.subplots(figsize=(7,7), dpi=100, nrows=1, ncols=2)
        self.fig1, self.ax1 = plt.subplots(figsize=(4,3), dpi=100, sharex=True, nrows=2, ncols=2)
        self.fig2, self.ax2 = plt.subplots(figsize=(6,6), dpi=100, sharex=True, nrows=2, ncols=2)
        self.ax[1].set_yticklabels([])
        self.fig.subplots_adjust( wspace=0.00, hspace=0.0,left=0.06, right=0.98, top=0.99, bottom=0.07)
        self.ax = self.ax.flatten()
        self.ax1 = self.ax1.flatten()
        self.ax2 = self.ax2.flatten()
        
        # Init images
        self.im = [None]*2
        self.im[0] = self.ax[0].imshow(im.histo_opt(self.m.temp[self.t0,:,:,self.itau],self.thres), cmap='gist_heat')
        self.im[1] = self.ax[1].imshow(im.histo_opt(self.m.vlos[self.t0,:,:,self.itau],self.thres), cmap='gist_gray')


        # Init fits
        self.im1_0 = [None]*self.ax1.size
        self.im1_1 = [None]*self.ax1.size

        for ii in range(4):
            self.im1_0[ii], = self.ax1[ii].plot(self.w[self.idx], self.o.dat[self.t0,self.yy, self.xx, self.idx, ii].squeeze(), 'ko', markersize=2)
            self.im1_1[ii], = self.ax1[ii].plot(self.w, self.s.dat[self.t0,self.yy, self.xx, :, ii].squeeze(), '-', color='orangered')

        self.ax1[1].set_ylim(-0.03,0.03); self.ax1[2].set_ylim(-0.03,0.03);self.ax1[3].set_ylim(-0.1,0.1)

        # Init model plots
        self.im2 = [None]*6
        self.im3 = [None]*4

        ix = [self.m.ltau[self.t0,self.yy, self.xx,self.itau],self.m.ltau[self.t0,self.yy, self.xx,self.itau]]
        
        self.im2[0], = self.ax2[0].plot(self.m.ltau[self.t0,self.yy, self.xx,:],self.m.temp[self.t0,self.yy, self.xx,:],'k-')
        self.im3[0], = self.ax2[0].plot(ix, [2000, 15000], color='orangered', linewidth=1)
        self.ax2[0].set_ylim(3500, 10000)

        self.im2[1], = self.ax2[1].plot(self.m.ltau[self.t0,self.yy, self.xx,:],self.m.vlos[self.t0,self.yy, self.xx,:], color='k')
        self.im2[2], = self.ax2[1].plot(self.m.ltau[self.t0,self.yy, self.xx,:],self.m.vturb[self.t0,self.yy, self.xx,:], color='dodgerblue')
        self.im3[1], = self.ax2[1].plot(ix, [-15, 15], color='orangered', linewidth=1)
        self.ax2[1].set_ylim(-7, 40)

        self.im2[3], = self.ax2[2].plot(self.m.ltau[self.t0,self.yy, self.xx,:],self.m.Bln[self.t0,self.yy, self.xx,:],'k-')
        self.im3[2], = self.ax2[2].plot(ix, [0, 4000], color='orangered', linewidth=1)
        self.ax2[2].set_ylim(-1500,1500)
        
        self.im2[4], = self.ax2[2].plot(self.m.ltau[self.t0,self.yy, self.xx,:],self.m.Bho[self.t0,self.yy, self.xx,:], color='k')
        self.im2[5], = self.ax2[3].plot(self.m.ltau[self.t0,self.yy, self.xx,:],self.m.azi[self.t0,self.yy, self.xx,:], color='dodgerblue')
        self.im3[3], = self.ax2[3].plot(ix, [0,180], color='orangered', linewidth=1)
        self.ax2[3].set_ylim(0,180)
        
        #start adding buttons and text
        rowcounter=1

    
        self.frame.pack()

        # radio buttons for display var
        profile_label = Tk.Label(self.frame, text='Panel 1:')
        profile_label.grid(row=rowcounter, column=0)
        self.p1v = Tk.StringVar()
        self.p1v.set("T")
        self.p2v = Tk.StringVar()
        self.p2v.set("V")
        self.pslider = Tk.DoubleVar()
        self.pslider.set(self.m.ltau[0,self.yy, self.xx,self.itau])
        
        i=1
        for option in self.names:
            r = Tk.Radiobutton(self.frame,text=option, variable=self.p1v, 
                               value=option[0],command=self.redrawIm1)
            r.grid(row=rowcounter,column=i,padx=0)
            i+=1
        rowcounter+=1

        profile_label1 = Tk.Label(self.frame, text='Panel 2:')
        profile_label1.grid(row=rowcounter, column=0)
        i=1
        for option in self.names:
            r = Tk.Radiobutton(self.frame,text=option, variable=self.p2v, 
                               value=option[0],command=self.redrawIm2)
            r.grid(row=rowcounter,column=i,padx=0)
            i+=1
        rowcounter+=1
        
        
        I0_label = Tk.Label(self.frame, text='Optical depth:').grid(row=rowcounter, column=0)
        I0_entry = Tk.Scale(self.frame,
                            orient='horizontal',
                            length=slider_length,
                            from_=self.m.ltau.min(), 
                            to=self.m.ltau.max(),
                            resolution=0.05,
                            variable=self.pslider,
                            command=self.updateDepth,
                            bigincrement=1.0)
        I0_entry.grid(row=rowcounter, column=1, columnspan=2)
        rowcounter+=1
        
        # quit button
        quit_button = Tk.Button(self.frame,text='Quit',command=self.quit)
        quit_button.grid(row=rowcounter, column=0, columnspan=2,sticky='ew')
        rowcounter+=1

        # Listen to clicks in x-Y
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.getXY)
        self.multi = MultiCursor(self.fig.canvas, (self.ax[0], self.ax[1]), color='white', lw=0.75,
                                 useblit=True, vertOn=True, horizOn=True)

        
        self.fig.show()
        self.fig1.set_tight_layout(True)
        self.fig1.show()
        self.fig2.set_tight_layout(True)
        self.fig2.show()

    def updateDepth(self, event):
        bla = float(self.pslider.get())
        self.itau = np.argmin(np.abs(self.m.ltau[self.t0,self.yy, self.xx,:]-bla))
        self.redrawIm1()
        self.redrawIm2()
        
        ix = [self.m.ltau[self.t0,self.yy, self.xx,self.itau],self.m.ltau[self.t0,self.yy, self.xx,self.itau]]

        for kk in range(4):
            self.im3[kk].set_xdata(ix)
        self.fig2.canvas.draw()

    def quit(self):
        plt.close("all")
        self.master.destroy()


    def getXY(self, event):
        if(event.inaxes == self.ax[0] or event.inaxes == self.ax[1]):
            self.xx = int(event.xdata+0.5)
            self.yy = int(event.ydata+0.5)
            print("(x,y) = ({0},{1})".format(self.xx, self.yy))
            self.reDrawPlots()
            #print(self.xx, self.yy)

    def redrawIm1(self):
        bla = str(self.p1v.get())
        if(bla == "T"):   self.v1 = self.vars[0]
        elif(bla == "V"): self.v1 = self.vars[1]
        elif(bla == "M"): self.v1 = self.vars[2]
        elif(bla == "B"): self.v1 = self.vars[3]
        elif(bla == "A"): self.v1 = self.vars[4]
        elif(bla == "I"): self.v1 = self.vars[5]
        elif(bla == "P"): self.v1 = self.vars[6]
        elif(bla == "B"): self.v1 = self.vars[7]

        ble = im.histo_opt(np.nan_to_num(self.v1[self.t0,:,:,self.itau]), self.thres)
        self.im[0].set_data(ble)
        self.im[0].set_clim(vmin=ble.min(),vmax=ble.max())
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        #self.fig.canvas.blit(self.ax[0].bbox)
        #stop()


    def redrawIm2(self):
        bla = str(self.p2v.get())
        if(bla == "T"):   self.v2 = self.vars[0]
        elif(bla == "V"): self.v2 = self.vars[1]
        elif(bla == "M"): self.v2 = self.vars[2]
        elif(bla == "L"): self.v2 = self.vars[3]
        elif(bla == "A"): self.v2 = self.vars[4]
        elif(bla == "H"): self.v2 = self.vars[5]
        elif(bla == "P"): self.v2 = self.vars[6]
        elif(bla == "B"): self.v2 = self.vars[7]

        ble = im.histo_opt(np.nan_to_num(self.v2[self.t0,:,:,self.itau]), self.thres)
        self.im[1].set_data(ble)
        self.im[1].set_clim(vmin=ble.min(),vmax=ble.max())
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

        
    def reDrawPlots(self):
        
        for ii in range(4):
            self.im1_0[ii].set_data(self.w[self.idx],self.o.dat[self.t0,self.yy, self.xx, self.idx, ii].squeeze())
            self.im1_1[ii].set_data(self.w,self.s.dat[self.t0,self.yy, self.xx, :, ii].squeeze())

            if(ii == 0):
                self.ax1[ii].relim()
                self.ax1[ii].autoscale_view()

            #self.fig1.canvas.draw()
            #self.fig1.canvas.blit(self.ax1[ii].bbox)
        self.fig1.canvas.draw()
        #self.fig1.canvas.flush_events()

            # Update model plots
        itau = self.m.ltau[self.t0,self.yy,self.xx,:]
        self.im2[0].set_data(itau,self.m.temp[self.t0,self.yy, self.xx, :])
        self.im2[1].set_data(itau,self.m.vlos[self.t0,self.yy, self.xx, :])
        self.im2[2].set_data(itau,self.m.vturb[self.t0,self.yy, self.xx, :])
        self.im2[3].set_data(itau,self.m.Bln[self.t0,self.yy, self.xx, :])
        self.im2[4].set_data(itau,self.m.Bho[self.t0,self.yy, self.xx, :])
        self.im2[5].set_data(itau,self.m.azi[self.t0,self.yy, self.xx, :])
        
        #self.fig2.canvas.draw()
        kk=0
        ix = [self.m.ltau[self.t0,self.yy, self.xx,self.itau],self.m.ltau[self.t0,self.yy, self.xx,self.itau]]

        for ii in self.ax2:
            self.im3[kk].set_xdata(ix)
            #self.fig2.canvas.blit(ii.bbox)
            kk+=1
        self.fig2.canvas.draw()
        #self.fig2.canvas.flush_events()
        

        
if __name__ == "__main__":
    root = Tk.Tk()
    app = iGui(root,'observed.nc','synthetic_cycle341.nc','atmosout_restart341.nc', t0=9)
#    app = iGui(root,'observed.nc','synthetic_cycle31.nc','atmosout_restart311.nc')

    root.mainloop()
    
