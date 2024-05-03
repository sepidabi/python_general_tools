import numpy as np
import os
import scipy.io.idl as idl
import astropy.units as u
import astropy.constants as const

class atlas:
    """
    Class to load (FTS) spectral atlas

    Methods:
        __init__()
        to(usys_to, perHz=True)
        get(w0, w1, cgs=False, si=False, nograv=False)

    Example:
        import atlas as S
        fts = S.atlas()
        wav, sp, cont = fts.get(6562.,6564., cgs=True, perHz=False)

    Author:
        Jaime de la Cruz Rodriguez (ISP/SU 2019)
    """
    def __init__(self):
        # Check dir where this class is stored
        this_dir, this_filename = os.path.split(__file__)
        DATA_PATH = os.path.join(this_dir, "../data/fts_disk_center_SI.idlsave")

        # Load data file
        fts = idl.readsav(DATA_PATH)
        self.cont = np.copy(fts["ftscnt_SI"])
        self.spec = np.copy(fts["ftsint_SI"])
        self.wave = np.copy(fts["ftswav"])
        self.usys = 'si_inu' # J/s/m^2/sr/Hz
        self.sunit = u.J / (u.s * u.m**2 * u.steradian * u.Hz)
        self.wunit = u.Angstrom


    def to(self, usys_to, perHz=True):
        usys_from = self.usys.lower()

        # Determine SI <-> cgs conversion
        if usys_to.lower() == 'cgs' and usys_from[:2] == 'si':
            conversion = u.J.to('erg') / (u.m.to('cm')**2)
            self.sunit *= u.m**2 / u.J * u.erg / u.cm**2
        elif usys_to.lower() == 'si' and usys_from[:3] == 'cgs':
            conversion = u.erg.to('J') / (u.cm.to('m')**2)
            self.sunit *= u.cm**2 / u.erg * u.J / u.m**2
        else:
            conversion = 1.

        # Apply I_lambda (per AA) <-> I_nu (per Hz) if need be
        lambda_to_nu = (self.wave*u.Angstrom.to('m'))**2 / const.c.value
        if (perHz == False and usys_from[-3:] != 'inu') or \
            (perHz == True and usys_from[-3:] == 'inu'):
            # no change to conversion factor
            if perHz == True:
                ext = '_inu'
            else:
                ext = '_ilambda'
        elif (perHz == False and usys_from[-3:] == 'inu'):
            conversion /= lambda_to_nu
            self.sunit *= u.Hz / u.Angstrom
            ext = '_ilambda'
        else:
            conversion *= lambda_to_nu
            self.sunit *= u.Angstrom / u.Hz
            ext = '_inu'

        # Apply conversion and update current unit system
        self.spec *= conversion
        self.cont *= conversion
        self.usys = usys_to + ext

    def get(self, w0, w1, cgs = False, si = False, nograv = False, perHz=True):
        idx = (np.where((self.wave >= w0) & (self.wave <= w1)))[0]
        
        if(cgs):
            self.to('cgs', perHz=perHz)
        elif(si):
            self.to('si', perHz=perHz)

        wave = np.copy(self.wave[idx[0]:idx[-1]])
        spec = np.copy(self.spec[idx[0]:idx[-1]])
        cont = np.copy(self.cont[idx[0]:idx[-1]])

        if(not nograv):
            wave *=  (1.0-633.0/const.c.value) # grav reddening

        # Normalize by the continuum if cgs=False and si=False (default)
        if (not cgs and not si):
            spec /= cont
            cont[:] = 1.0
            
        return wave, spec, cont
