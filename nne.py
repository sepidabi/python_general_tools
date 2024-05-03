import nicole as nic
import sparsetools as sp
import numpy as np


m = sp.model('/srv/scratch/crobu/sst/2013.07.15/inv8542_111_114/inv/modelout_7_z.nc')

a=m.nne[0,:,:,:].transpose()

a.tofile('/srv/scratch/crobu/sst/2013.07.15/inv8542_111_114/inv/ne')

