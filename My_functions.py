import numpy as np
import lptools as lp
from sklearn.cluster import KMeans
import matplotlib.pyplot as ply


def read_polcube(fil0,fil1):

    nx, ny, dum, ns, dtype, ndim = lp.lphead(fil0)
    nw, nt, dum, ns, dtype, ndim = lp.lphead(fil1)

    io = np.memmap(fil0, mode='r', shape=(nt,ns,nw,ny,nx), dtype=dtype,
               offset=512)
    return io
    
def read_unpolcube(fil0,fil1):

    nx, ny, dum, dtype, ndim = lp.lphead(fil0)
    nw, nt, dum, dtype, ndim = lp.lphead(fil1)

    io = np.memmap(fil0, mode='r', shape=(nt,nw,ny,nx), dtype=dtype,
               offset=512)
    return io

def read_lp_cube(file):
    nx, ny, nt, dum, dtype, ndim = lp.lphead(file)
    io = np.memmap(file, mode='r', shape=(nt,ny,nx), dtype=dtype,
               offset=512)
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



def derivative(x,y):

    n= x.shape[0]
    yp= np.zeros(n)
    
    dx  = x[1:n-1] -x[0:n-2]
    dx1 = x[2:] - x[1:n-1]

    der = (y[1:n-1] - y[0:n-2]) / dx
    der1 = (y[2:] - y[1:n-1]) / dx1

    idx= np.where(der*der1 > 0.0)
    if (((der*der1 > 0.0).nonzero() >0) == True):
        lambdad = (1. + dx1[idx] / (dx1[idx] + dx[idx])) / 3.
        yp[np.asarray(idx)+1] = (der[idx] * der1[idx]) / (lambdad * der1[idx] + (1. - lambdad) * der[idx]) 
    yp[0] = der[0]  
    yp[n-1] = der1[der1.shape[0]-1]

    return yp



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
