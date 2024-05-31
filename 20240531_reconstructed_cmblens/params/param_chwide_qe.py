"""
Parameter file for S4 chwide simulations 

For now we focus on the QE only, 
but in principle this param file can easly be extended to the MAP estimator.
"""

import os
from os.path import join as opj
import numpy as np
import healpy as hp

import plancklens

from plancklens import utils, qresp, qest, qecl
from plancklens.qcinv import cd_solve
from plancklens.filt import filt_cinv, filt_util
from plancklens.utils import cli
from plancklens import nhl
from plancklens.n1 import n1

import os
import sys
module_path = '/global/homes/l/llegrand/litebirdXS4/20240531_reconstructed_cmblens/sims/'
if module_path not in sys.path:
    sys.path.append(module_path)
import chwide_sims


suffix = 'chwide' # descriptor to distinguish this parfile from others...
TEMP = opj(os.environ['SCRATCH'], 'lenscarfrecs', suffix)
    
lmax_cmb_len = 4096 
lmax_cmb_unl = 4096 + 1024

# TODO: Check if that value is correct 
lmin_cmb_transf = 1 

#----------------- Lensing reconstruction info
lmax_ivf, mmax_ivf, beam, nlev_t, nlev_p = (lmax_cmb_len, lmax_cmb_len, 2.1, 1., 1.3)

lmin_tlm, lmin_elm, lmin_blm = (lmin_cmb_transf, lmin_cmb_transf, lmin_cmb_transf) # The fiducial transfer functions are set to zero below these lmins

lmax_qlm, mmax_qlm = (lmax_cmb_len, lmax_cmb_len) # Lensing map is reconstructed down to this lmax and mmax
# NB: the QEs from plancklens does not support mmax != lmax, but the MAP pipeline does                                                                                                                                                                                                                         the MAP pipeline does


# Multigrid chain descriptor
# The hard coded number nside 2048 here is irrelevant for diagonal preconditioner
chain_descrs = lambda lmax_sol, cg_tol : [[0, ["diag_cl"], lmax_sol, nside, np.inf, cg_tol, cd_solve.tr_cg, cd_solve.cache_mem()]]
libdir_iterators = lambda qe_key, simidx, version, tol: opj(TEMP,'%s_sim%04d'%(qe_key, simidx) + version + '_tol' *(tol != 5))
#------------------

# Fiducial CMB spectra for QE and iterative reconstructions
# (here we use very lightly suboptimal lensed spectra QE weights)
cls_path = opj(os.path.dirname(plancklens.__file__), 'data', 'cls')
cls_unl = utils.camb_clfile(opj(cls_path, 'FFP10_wdipole_lenspotentialCls.dat'))
cls_len = utils.camb_clfile(opj(cls_path, 'FFP10_wdipole_lensedCls.dat'))
cls_grad = utils.camb_clfile(opj(cls_path, 'FFP10_wdipole_gradlensedCls.dat'))
cls_weight = cls_grad 

# Fiducial model of the transfer function
transf_tlm   =  hp.gauss_beam(beam/180 / 60 * np.pi, lmax=lmax_ivf) * (np.arange(lmax_ivf + 1) >= lmin_tlm)
transf_elm   =  hp.gauss_beam(beam/180 / 60 * np.pi, lmax=lmax_ivf) * (np.arange(lmax_ivf + 1) >= lmin_elm)
transf_blm   =  hp.gauss_beam(beam/180 / 60 * np.pi, lmax=lmax_ivf) * (np.arange(lmax_ivf + 1) >= lmin_blm)
transf_d = {'t':transf_tlm, 'e':transf_elm, 'b':transf_blm}

# Isotropic approximation to the filtering (used eg for response calculations)
ftl = cli(cls_len['tt'][:lmax_ivf + 1] + (nlev_t / 180 / 60 * np.pi) ** 2 * cli(transf_tlm ** 2)) * (transf_tlm > 0)
fel = cli(cls_len['ee'][:lmax_ivf + 1] + (nlev_p / 180 / 60 * np.pi) ** 2 * cli(transf_elm ** 2)) * (transf_elm > 0)
fbl = cli(cls_len['bb'][:lmax_ivf + 1] + (nlev_p / 180 / 60 * np.pi) ** 2 * cli(transf_blm ** 2)) * (transf_blm > 0)

# Same using unlensed spectra (used for unlensed response used to initiate the MAP curvature matrix)
ftl_unl = cli(cls_unl['tt'][:lmax_ivf + 1] + (nlev_t / 180 / 60 * np.pi) ** 2 * cli(transf_tlm ** 2)) * (transf_tlm > 0)
fel_unl = cli(cls_unl['ee'][:lmax_ivf + 1] + (nlev_p / 180 / 60 * np.pi) ** 2 * cli(transf_elm ** 2)) * (transf_elm > 0)
fbl_unl = cli(cls_unl['bb'][:lmax_ivf + 1] + (nlev_p / 180 / 60 * np.pi) ** 2 * cli(transf_blm ** 2)) * (transf_blm > 0)


# -------------------------
# ---- Input simulation libraries

sims = chwide_sims.chwide_nilc()
nside = sims.nside

# ------------------------- Masks

# List of paths to masks that will be multiplied together to give the total mask
# Here we add together the E and B masks 
# In practice E maps have a different mask as the B one, not sure if we can have a different mask for each

mask_path = '/global/cfs/cdirs/cmbs4xlb/v1/component_separated/chwide/masks_common'
masks = [opj(mask_path, 'dust-mask-15pc_3pc-apo_NSIDE2048.fits'), opj(mask_path, 'chwide_clip0p3relhits_3degC2apo_NSIDE2048.fits')]


# ------------------------- Filtering

# List of the inverse noise pixel variance maps, all will be multiplied together
ninv_t = [np.array([hp.nside2pixarea(nside, degrees=True) * 60 ** 2 / nlev_t ** 2])] + masks
cinv_t = filt_cinv.cinv_t(opj(TEMP, 'cinv_t'), lmax_ivf, nside, cls_len, transf_tlm, ninv_t,
                        marge_monopole=True, marge_dipole=True, marge_maps=[])

ninv_p = [[np.array([hp.nside2pixarea(nside, degrees=True) * 60 ** 2 / nlev_p ** 2])] + masks]
cinv_p = filt_cinv.cinv_p(opj(TEMP, 'cinv_p'), lmax_ivf, nside, cls_len, transf_elm, ninv_p,
            chain_descr=chain_descrs(lmax_ivf, 1e-5), transf_blm=transf_blm, marge_qmaps=(), marge_umaps=())

# THE CLS WEIGHTS ARE THE GRAD CLS
ivfs_raw    = filt_cinv.library_cinv_sepTP(opj(TEMP, 'ivfs'), sims, cinv_t, cinv_p, cls_weight)
ftl_rs = np.ones(lmax_ivf + 1, dtype=float) * (np.arange(lmax_ivf + 1) >= lmin_tlm)
fel_rs = np.ones(lmax_ivf + 1, dtype=float) * (np.arange(lmax_ivf + 1) >= lmin_elm)
fbl_rs = np.ones(lmax_ivf + 1, dtype=float) * (np.arange(lmax_ivf + 1) >= lmin_blm)
ivfs   = filt_util.library_ftl(ivfs_raw, lmax_ivf, ftl_rs, fel_rs, fbl_rs)

# ---- QE libraries from plancklens to calculate unnormalized QE (qlms) and their spectra (qcls)
mc_sims_bias = np.arange(60, dtype=int)
mc_sims_var  = np.arange(60, 300, dtype=int)
qlms_dd = qest.library_sepTP(opj(TEMP, 'qlms_dd'), ivfs, ivfs,  cls_weight['te'], nside, lmax_qlm=lmax_qlm)
qcls_dd = qecl.library(opj(TEMP, 'qcls_dd'), qlms_dd, qlms_dd, mc_sims_bias)


# -------------------------
# This following block is only necessary if a full, Planck-like QE lensing power spectrum analysis is desired
# This uses 'ds' and 'ss' QE's, crossing data with sims and sims with other sims.

# This remaps idx -> idx + 1 by blocks of 60 up to 300. This is used to remap the sim indices for the 'MCN0' debiasing term in the QE spectrum
ss_dict = { k : v for k, v in zip( np.concatenate( [ range(i*60, (i+1)*60) for i in range(0,5) ] ),
                                   np.concatenate( [ np.roll( range(i*60, (i+1)*60), -1 ) for i in range(0,5) ] ) ) }
ds_dict = { k : -1 for k in range(300)} # This remap all sim. indices to the data maps to build QEs with always the data in one leg

ivfs_d = filt_util.library_shuffle(ivfs, ds_dict)
ivfs_s = filt_util.library_shuffle(ivfs, ss_dict)

qlms_ds = qest.library_sepTP(opj(TEMP, 'qlms_ds'), ivfs, ivfs_d, cls_weight['te'], nside, lmax_qlm=lmax_qlm)
qlms_ss = qest.library_sepTP(opj(TEMP, 'qlms_ss'), ivfs, ivfs_s, cls_weight['te'], nside, lmax_qlm=lmax_qlm)

qcls_ds = qecl.library(opj(TEMP, 'qcls_ds'), qlms_ds, qlms_ds, np.array([]))  # for QE RDN0 calculations
qcls_ss = qecl.library(opj(TEMP, 'qcls_ss'), qlms_ss, qlms_ss, np.array([]))  # for QE RDN0 / MCN0 calculations

# -------------------------
# QE response and N0 and N1 biases  

libdir_nhl_dd = os.path.join(TEMP, 'nhl_dd')
nhl_dd = nhl.nhl_lib_simple(libdir_nhl_dd, ivfs, cls_weight, lmax_qlm)

libdir_n1_dd = os.path.join(TEMP, 'n1_ffp10')
n1_dd = n1.library_n1(libdir_n1_dd,cls_len['tt'],cls_len['te'],cls_len['ee'])
# The cls_len are the one in the map. The cls_weights should be given in the arguments of get_n1

libdir_resp_dd = os.path.join(TEMP, 'qresp')
qresp_dd = qresp.resp_lib_simple(libdir_resp_dd, lmax_ivf, cls_weight, cls_grad,
                                 {'t': ivfs.get_ftl(), 'e':ivfs.get_fel(), 'b':ivfs.get_fbl()}, lmax_qlm)
