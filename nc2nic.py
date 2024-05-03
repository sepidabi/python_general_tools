import nicole as nic
import sparsetools as sp
import numpy as np


m = sp.model('/srv/scratch/crobu/sst/2013.07.15/inv8542_111_114/inv/modelout_1_z.nc')

m0 = nic.model('/srv/scratch/crobu/sst/2013.07.15/inv8542_111_114/inv/modelout_1.nic')

m0.z[:,:,:]=m.z[0,:,:,:]


m0.write('/srv/scratch/crobu/sst/2013.07.15/inv8542_111_114/inv/modelout_1_z.nic')


