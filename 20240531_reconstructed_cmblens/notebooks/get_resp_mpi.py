"""
Simple script to get the effective CMB lensing reponse with MPI 
"""
import numpy as np 
from plancklens import utils
import healpy as hp
from plancklens.helpers import mpi
from plancklens.helpers import cachers
from os.path import join as opj
import sys 
import os

print('Running with rank %s of %s' % (mpi.rank, mpi.size))


# Import param file
module_path = os.path.abspath(os.path.join('../params'))
if module_path not in sys.path:
    sys.path.append(module_path)

import param_chwide_qe_lminB200_v2 as param_file

# Import sims library 

module_path = os.path.abspath(os.path.join('../../20231018_signal_cmb'))
if module_path not in sys.path:
    sys.path.append(module_path)
    
import unlensed_ffp10

cmb_unl= unlensed_ffp10.cmb_unl_ffp10()

k = 'p_p'  # Planck lensing QE key

savedir = '/global/cfs/cdirs/cmbs4xlb/v1/lensingrec'
tag = 'chwide_qe_v2'

fn_cache_cpp = opj(savedir, tag, 'cpp_cache')

cacher = cachers.cacher_npy(fn_cache_cpp)

fn_ckk_pxin = lambda simidx: f'ckk_pxin_{simidx:04d}'
fn_ckk_in = lambda simidx: 'ckk_in_%04d'%  (simidx)


nsims = 500
mc_sims = np.arange(nsims)
fsky = np.mean(param_file.ivfs.get_fmask())


def get_resp_eff(k, idx):
    if not cacher.is_cached(fn_ckk_in(idx)) or not cacher.is_cached(fn_ckk_pxin(idx)) or recache:
        print(f'rank {mpi.rank} computing sim {idx}')
        # resp = param_file.qresp_dd.get_response(k, 'p')
        qlms_dd = param_file.qlms_dd
        lmax = param_file.qresp_dd.lmax_qlm
        cpp_qexin = np.zeros(lmax+1)

        # for idx in mc_sims:
            # print(idx)
        plm_qe = qlms_dd.get_sim_qlm(k, idx) 
        plm_in = utils.alm_copy(cmb_unl.get_sim_plm(idx), lmax=hp.Alm.getlmax(np.size(plm_qe)))
        cpp_qexin = hp.alm2cl(plm_qe, plm_in) / fsky
        cpp_in  = hp.alm2cl(plm_in)
        cacher.cache(fn_ckk_in(idx), cpp_in)
        cacher.cache(fn_ckk_pxin(idx), cpp_qexin)
        # return cpp_qexin

recache = False

for idx in mc_sims[mpi.rank::mpi.size]:
    get_resp_eff(k, idx)


# # Save the response 
# fn_resp = os.path.join(savedir, tag, 'resp_qe.txt')
# if not os.path.exists(fn_resp):
#     resp_arr = np.array([np.arange(resp.shape[0]), resp])
#     np.savetxt(fn_resp, resp_arr.T, header='Fiducial QE response\n L, Resp')

# recache = True 

fn_resp = os.path.join(savedir, tag, 'resp_qe_effective.txt')
if mpi.rank == 0 and not os.path.exists(fn_resp) or recache:
    resp_eff = np.zeros(param_file.lmax_qlm + 1)
    for idx in mc_sims:
        cpp_qexin = cacher.load(fn_ckk_pxin(idx))
        cpp_in = cacher.load(fn_ckk_in(idx))
        resp_eff += cpp_qexin * utils.cli(cpp_in)
    resp_eff /= len(mc_sims)
    resp_arr = np.array([np.arange(resp_eff.shape[0]), resp_eff])
    np.savetxt(fn_resp, resp_arr.T, header='Effective QE response\n L, Resp')

