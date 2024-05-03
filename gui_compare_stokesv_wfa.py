import matplotlib
matplotlib.use("TkAgg")
import Tkinter as Tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
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
    def __init__(self, datamicrojet, dstokesi, wav):
        slider_length = 250
        self.mic = idlsave.read(datamicrojet)
        self.di = idlsave.read(dstokesi)
        self.wav = idlsave.read(wav)

        self.x = self.mic.xx - self.mic.coord[0]
        self.y = self.mic.yy - self.mic.coord[1]
        #self.t = self.mic.tt

        self.wavcal = (self.wav.spect_pos - self.wav.spect_pos[10]) #mA

        self.t = input("Frame: ")
        #self.t = float(self.t)
        
        #self.nframes = np.arange(self.mic.frame_f - self.mic.frame_i + 1) #+ self.mic.frame_i
        
        self.master = root
        self.master.wm_title("Control window")
        self.frame = Tk.Frame(self.master)

        plt.ion()

        #self.fig = plt.figure(figsize=(15,7), dpi = 100)
        #self.ax01 = self.fig.add_subplot(121)
        #self.ax02 = self.fig.add_subplot(122)
        self.fig, self.ax = plt.subplots(figsize=(15,7), nrows = 1, ncols = 2, dpi = 100)
        self.fig1, self.ax1 = plt.subplots(figsize=(8,5), nrows = 1, ncols = 1, dpi = 100)
        self.fig2, self.ax2 = plt.subplots(figsize=(8,5), nrows = 1, ncols = 1, dpi = 100)
        self.ax3 = self.ax2.twinx()
        self.fig.subplots_adjust( wspace=0.00, hspace=0.0,left=0.06, right=0.98, top=0.99, bottom=0.07)
        self.ax = self.ax.flatten()

        self.im01 = [None]
        self.im01 = self.ax[0].imshow(im.histo_opt(self.mic.profileca[0, self.t, 6, :, :]), cmap = 'gist_gray', origin = 'lower')
        self.im02 = [None]
        self.im02 = self.ax[1].imshow(im.histo_opt(self.mic.profileca[3, self.t, 16, :, :]), cmap = 'gist_gray', origin = 'lower')
                                       
        self.im1, = [None]
        self.im1, = self.ax1.plot(self.wavcal, self.mic.profileca[0, self.t, :, self.y, self.x], 'k-')
        
        self.im2, = [None]
        self.im2, = self.ax2.plot(self.wavcal, self.mic.profileca[3, self.t, :, self.y, self.x].squeeze(), 'k-')
        self.im3, = [None]
        self.im3, = self.ax3.plot(self.wavcal, self.di.di[:, self.y, self.x].squeeze(), 'm-')
        
        rowcounter=1

        self.frame.pack()

        #profile_label = Tk.Label(self.frame, text='Panel 1:')
        #profile_label.grid(row=rowcounter, column=0)
        #self.p1v = Tk.StringVar()
        #self.p2v = Tk.StringVar()
        #self.pslider = Tk.DoubleVar()
        #self.pslider.set(self.t)
        
        #I0_label = Tk.Label(self.frame, text='Frame:').grid(row=rowcounter, column=0)
        #I0_entry = Tk.Scale(self.frame,orient='horizontal', length=slider_length, from_=self.nframes.min(), to=self.nframes.max(),resolution=1.,variable=self.pslider,command=self.updateFrame)
        #I0_entry.grid(row=rowcounter, column=1, columnspan=2)
        #rowcounter+=1


        left_button = Tk.Button(self.frame,text='<',width = 10, command = self.leftframe)
        #left_button.pack(side=left)
        left_button.grid(row=2, column=0, padx = 5, pady = 5)

        #self.ts = str(self.t)
        self.ts = Tk.StringVar() #str(self.t)
        label = Tk.Label(self.frame, textvariable = self.ts, font  = ("Helvetica", 12))
        label.grid(row = 2, column = 1)
        
        right_button = Tk.Button(self.frame,text='>',width = 10, command = self.rightframe)
        #right_button.pack(side=right)
        right_button.grid(row=2, column=3, padx = 5, pady = 5)
        
        quit_button = Tk.Button(self.frame,text='Quit',command=self.quit, width = 10)
        quit_button.grid(row=4, column=0, columnspan = 4, sticky = 'ew', padx = 5, pady = 5)
        #rowcounter+=1

        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.getXY)
        self.multi = MultiCursor(self.fig.canvas, (self.ax[0], self.ax[1]), color='w', lw=0.75,useblit=True, vertOn=True, horizOn=True)

        self.fig.show()
        self.fig1.set_tight_layout(True)
        self.fig1.show()
        self.fig2.set_tight_layout(True)

        
    #def updateFrame(self, event):
        #bla = float(self.pslider.get())
        #self.t = np.argmin(np.abs(bla))
        #self.reDrawPlots()
        
        #ix = [self.t,self.t]

        #self.im2.set_xdata(ix)
        #self.fig2.canvas.draw()
    

    def rightframe(self):
        self.t +=1
        print(self.t)
        self.reDrawPlots()
        self.reShowImages()
        self.ts.set(self.t)

    def leftframe(self):
        self.t -=1
        print(self.t)
        self.reDrawPlots()
        self.reShowImages()
        self.ts.set(self.t)
        
    def quit(self):
        plt.close("all")
        self.master.destroy()

    def getXY(self, event):
        if(event.inaxes == self.ax[0] or event.inaxes == self.ax[1]):
            self.x = int(event.xdata)
            self.y = int(event.ydata)
            print("(x,y) = ({0},{1})".format(self.x, self.y))
            self.reDrawPlots()

    def reShowImages(self):
        self.fig.canvas.draw()
        self.im01.set_data(self.mic.profileca[0, self.t, 6, :, :])
        self.im02.set_data(self.mic.profileca[3, self.t, 16, :, :])
        
        
    def reDrawPlots(self):

        iframe = self.t
        self.fig1.canvas.draw()
        print(self.t)
        self.im1.set_data(self.wavcal, self.mic.profileca[0, self.t, :, self.y, self.x])
        self.fig1.canvas.flush_events()
        self.fig2.canvas.draw()
        self.im2.set_data(self.wavcal, self.mic.profileca[3, self.t, :, self.y, self.x])
        self.im3.set_data(self.wavcal, self.di.di[:, self.y, self.x])
        self.fig2.canvas.flush_events()
        #self.im1.set_data(self.wavcal, self.mic.profileca[0,self.t-self.mic.frame_i,:,self.y, self.x].squeeze())
        #self.im2.set_data(self.wavcal, self.mic.profileca[3,self.t-self.mic.frame_i,:,self.y, self.x].squeeze(), 'k-')
        #self.im3.set_data(self.wavcal, di.di[:,self.y, self.x].squeeze(), 'm-')
            
        #if(ii == 0):
        #    self.ax1[ii].relim()
        #    self.ax1[ii].autoscale_view()

            #self.fig1.canvas.draw()
            #self.fig1.canvas.blit(self.ax1[ii].bbox)
        self.fig1.canvas.draw()
        self.fig1.canvas.flush_events()

        
if __name__ == "__main__":
    root = Tk.Tk()
    app = iGui('microjet_1_170515.sav','dI_microjet_1_170403.sav','/srv/scratch/seste/20130722_chromo/spectfile.8542.idlsave')

    root.mainloop()
    
