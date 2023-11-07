import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from .config import CMBDIR
import os 


class CMBLensed:
    """
    class for getting lensed CMB simulation
    """
    def __init__(self,nside=512)
        self.nside = nside

    def TEB(self,idx):
        fname =  os.path.join(CMBDIR,f'lcdm_teb_{idx:04d}.npy')
        return np.load(fname)

    def TQU(self,idx):
        return hp.alm2map(self.TEB(idx),self.nside)

    def Kappa(self,idx):
        fname = os.path.join(CMBDIR,f'lcdm_k_{idx:04d}.npy')
        return np.load(fname)

    def Phi(self,idx):
        raise NotImplementedError