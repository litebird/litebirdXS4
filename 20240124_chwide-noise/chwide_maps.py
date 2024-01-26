import numpy as np
UNSEEN = -1.6375e+30

try:
    from healpy import UNSEEN as hp_UNSEEN
    if UNSEEN != hp_UNSEEN:
        print("using %.5e for UNSEEN pixels" % UNSEEN)
except ImportError:
    UNSEEN = -1.6375e+30
    hp_UNSEEN = None
    print("using %.5e for UNSEEN pixels"%UNSEEN)

tht_min, tht_max = (64 / 180 * np.pi, 155. / 180 * np.pi)
pix_min, npix = 14135296, 33839732
path = '/global/cfs/cdirs/cmbs4xlb/v1/noise/chwide/tqu_%03d_%04d.npy'

class noise:
    @staticmethod
    def get_sim_tmap(freq, idx):
        m = np.full( (12 * 2048 ** 2), UNSEEN)
        m[pix_min: pix_min + npix] = np.load(path%(freq, idx))[0]
        return m
    @staticmethod
    def get_sim_qumap(freq, idx):
        m = np.full((2, 12 * 2048 ** 2), UNSEEN)
        m[:, pix_min:pix_min + npix] = np.load(path%(freq, idx))[1:]
        return m
    @staticmethod
    def get_sim_tqumap(freq, idx):
        m = np.full((3, 12 * 2048 ** 2), UNSEEN)
        m[:, pix_min:pix_min + npix] = np.load(path%(freq, idx))
        return m

