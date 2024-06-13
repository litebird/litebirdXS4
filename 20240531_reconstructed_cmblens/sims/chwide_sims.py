r"""S4 chwide  simulation libraries.

    Note:
        These simulations are located on NERSC systems.

    Note:
        Units of the maps stored at NERSC are :math:`K` but this module returns maps in :math:`\mu K`

"""

import os
from os.path import join as opj
import healpy as hp
import numpy as np

from plancklens import utils

has_key = lambda key : key in os.environ.keys()

#assert has_key('NERSC_HOST'), "Maps are stored on NERSC"

class chwide_nilc:
    r""" 
     Simulated chile wide polarization maps 

    """
    def __init__(self):
        chwide_dir = '/global/cfs/cdirs/cmbs4xlb/v1/component_separated/chwide'

        self.emaps = opj(chwide_dir, 'nilc_Emaps/fits', 'NILC_CMB-S4_CHWIDE-E%s_NSIDE2048_fwhm2.1_CHLAT-only_medium_cos-NSIDE2048-lmax4096_mc0%02d.fits' )
        
        self.bmaps = opj(chwide_dir, 'nilc_Bmaps/fits', 'NILC_CMB-S4_CHWIDE-B%s_NSIDE2048_fwhm2.1_CHLAT-only_medium_cos-NSIDE2048-lmax4096_galmask_mc0%02d.fits')

        self.lmax = 4096
        self.nside = 2048

    def hashdict(self):
        return {'emaps':self.emaps, 'bmaps':self.bmaps}

    def get_sim_elm(self, idx, map_type='map'):
        # assert map_type in ['map', 'fg', 'noise']
        emap =  hp.read_map(self.emaps % (map_type, idx))
        return hp.map2alm(emap, lmax=self.lmax)
    
    def get_sim_blm(self, idx, map_type='map'):
        # assert map_type in ['map', 'fg', 'noise']
        bmap =  hp.read_map(self.bmaps % (map_type, idx))
        return hp.map2alm(bmap, lmax=self.lmax)


    def get_sim_pmap(self,idx):
        """Returns polarization healpy maps for a simulation

            Args:
                idx: simulation index

            Returns:
                Q and U healpy maps

        """
        elm = self.get_sim_elm(idx)
        blm = self.get_sim_blm(idx)

        Q,U = hp.alm2map_spin([elm,blm], self.nside, 2, hp.Alm.getlmax(elm.size))
        del elm,blm
        return Q, U