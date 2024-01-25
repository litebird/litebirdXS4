import numpy as np
import healpy as hp

tht_min, tht_max = (64 / 180 * np.pi, 155. / 180 * np.pi)
pix_min, npix = 14135296, 33839732
path = '' #FIXME

class noise:
    @staticmethod
    def get_sim_tmap(freq, idx):
        m = np.full( (12 * 2048 ** 2), hp.UNSEEN)
        m[pix_min: pix_min + npix] = np.load(path%(freq, idx))[0]
        return m
    @staticmethod
    def get_sim_qumap(freq, idx):
        m = np.full((2, 12 * 2048 ** 2), hp.UNSEEN)
        m[:, pix_min:pix_min + npix] = np.load(path%(freq, idx))[1:]
        return m
    @staticmethod
    def get_sim_tqumap(freq, idx):
        m = np.full((3, 12 * 2048 ** 2), hp.UNSEEN)
        m[:, pix_min:pix_min + npix] = np.load(path%(freq, idx))
        return m

