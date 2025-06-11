import os
from os.path import join as opj
import numpy as np
import healpy as hp
from lenspyx import utils_hp

aberration_lbv_ffp10 = (264. * (np.pi / 180), 48.26 * (np.pi / 180), 0.001234)


class cmb_unl_ffp10:
    """FFP10 input sim libraries, unlensed alms.


        Note: unlensed blm are hdu=3 of these arrays, but are identically zero

    """
    def __init__(self, aberration:tuple[float, float, float]or None=None, verbose=True):
        if aberration is None:
            aberration = aberration_lbv_ffp10

        self.lmax_unl = 5120
        self.lmax_unl_tensor = 2048
        # aberration: we must add the difference to the FFP10 aberration
        l, b, v = aberration
        l_ffp10, b_ffp10, v_ffp10 = aberration_lbv_ffp10

        # \phi_{10} = - \sqrt{4\pi/3} n_z
        # \phi_{11} = + \sqrt{4\pi / 3} \frac{(n_x - i n_y)}{\sqrt{2}}
        vlm = np.array([0., np.cos(b), - np.exp(-1j * l) * np.sin(b) / np.sqrt(2.)])  # LM = 00, 10 and 11
        vlm_ffp10 = np.array([0., np.cos(b_ffp10), - np.exp(-1j * l_ffp10) * np.sin(b_ffp10) / np.sqrt(2.)])
        vlm       *= (-v * np.sqrt(4 * np.pi / 3))
        vlm_ffp10 *= (-v_ffp10 * np.sqrt(4 * np.pi / 3))
        self.delta_vlm = vlm - vlm_ffp10
        self.vlm = vlm
        if verbose:
            print("Input aberration power %.3e"%(utils_hp.alm2cl(vlm, vlm, 1, 1, 1)[1]))
    @staticmethod
    def hashdict():
        return {'sim_lib': 'ffp10 unlensed scalar cmb inputs'}

    @staticmethod
    def get_sim_tlm(idx):
        """
            Args:
                idx: simulation index

            Returns:
                unlensed temperature simulation healpy alm array

        """
        return 1e6 * hp.read_alm(opj(os.environ["CFS"],'cmb/data/generic/cmb/ffp10/mc/scalar/ffp10_unlensed_scl_cmb_000_tebplm_mc_%04d.fits'% idx), hdu=1)

    @staticmethod
    def get_sim_tensor_tlm(idx, r=0.01):
        """
            Args:
                idx: simulation index
                r: scalar to tensor ratio

            Returns:
                unlensed temperature tensor mode simulation healpy alm array for input r value

        """
        fn = opj(os.environ["CFS"],'cmb/data/generic/cmb/ffp10/mc/tensor/ffp10_ten_cmb_000_alm_mc_%04d.fits'%idx)
        return (1e6 * np.sqrt(r / 0.01)) * hp.read_alm(fn, hdu=1)

    @staticmethod
    def get_sim_tensor_elm(idx, r=0.01):
        """
            Args:
                idx: simulation index
                r: scalar to tensor ratio

            Returns:
                unlensed E-mode tensor mode simulation healpy alm array for input r value

        """
        fn = opj(os.environ["CFS"],'cmb/data/generic/cmb/ffp10/mc/tensor/ffp10_ten_cmb_000_alm_mc_%04d.fits'%idx)
        return (1e6 * np.sqrt(r / 0.01)) * hp.read_alm(fn, hdu=2)

    @staticmethod
    def get_sim_tensor_blm(idx, r=0.01):
        """
            Args:
                idx: simulation index
                r: scalar to tensor ratio

            Returns:
                unlensed B-mode tensor mode simulation healpy alm array for input r value

        """
        fn = opj(os.environ["CFS"],'cmb/data/generic/cmb/ffp10/mc/tensor/ffp10_ten_cmb_000_alm_mc_%04d.fits'%idx)
        return (1e6 * np.sqrt(r / 0.01)) * hp.read_alm(fn, hdu=3)

    @staticmethod
    def get_sim_elm(idx):
        """
            Args:
                idx: simulation index

            Returns:
                unlensed E-polarization simulation healpy alm array

        """
        return 1e6 * hp.read_alm(opj(os.environ["CFS"],'cmb/data/generic/cmb/ffp10/mc/scalar/ffp10_unlensed_scl_cmb_000_tebplm_mc_%04d.fits'% idx), hdu=2)


    def get_sim_plm(self, idx):
        r"""
            Args:
                idx: simulation index

            Returns:
               lensing potential :math:`\phi_{LM}` simulation healpy alm array

        """
        plm =  hp.read_alm(opj(os.environ["CFS"],'cmb/data/generic/cmb/ffp10/mc/scalar/ffp10_unlensed_scl_cmb_000_tebplm_mc_%04d.fits'% idx), hdu=4)
        lmax = hp.Alm.getlmax(plm.size)
        assert lmax == self.lmax_unl, (lmax, self.lmax_unl)
        plm[utils_hp.Alm.getidx(lmax, 1, 0)] += self.delta_vlm[1] # LM=10 aberration
        plm[utils_hp.Alm.getidx(lmax, 1, 1)] += self.delta_vlm[2] # LM = 11
        return plm
    def get_sim_dlm(self, idx):
        dlm = self.get_sim_plm(idx)
        lmax = hp.Alm.getlmax(dlm.size)
        assert lmax == self.lmax_unl, (lmax, self.lmax_unl)
        p2d = np.sqrt(np.arange(lmax + 1, dtype=float) * np.arange(1,lmax + 2, dtype=float))
        hp.almxfl(dlm, p2d, inplace=True)
        return dlm