import nicole as nic
import sparsetools as sp
import numpy as np



m0 = nic.model('/srv/scratch/crobu/sst/2013.07.15/inv8542_111_114/inv/modelout_1.nic')
m1 = sp.model(nx=m0.nx, ny=m0.ny, ndep=m0.nz, nt=1)

m1.temp[0,:,:,:] = m0.t
m1.pgas[0,:,:,:] = m0.gas_p
m1.ltau[0,:,:,:] = m0.tau
m1.vturb[0,:,:,:] = m0.v_mic
m1.vlos[0,:,:,:] = m0.vlos

m1.write('/srv/scratch/crobu/sst/2013.07.15/inv8542_111_114/inv/modelout_1.nc')

