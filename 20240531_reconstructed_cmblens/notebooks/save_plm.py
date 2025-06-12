
import numpy as np
import healpy as hp
import os
from plancklens import utils
import sys 
from plancklens.helpers import mpi

# Import param file
module_path = os.path.abspath(os.path.join('../params'))
if module_path not in sys.path:
    sys.path.append(module_path)

import param_chwide_qe_lminB200_v2 as param_file


""" Scripts to load and save the phi previously estimated
"""

qlms_dd = param_file.qlms_dd
k = 'p_p'
resp = param_file.qresp_dd.get_response('p_p', 'p')



nsims = 500
savedir = '/global/cfs/cdirs/cmbs4xlb/v1/lensingrec'
tag = 'chwide_qe_v2'
if not os.path.exists(os.path.join(savedir, tag)):
    os.makedirs(os.path.join(savedir, tag))

L, resp_eff = np.loadtxt(os.path.join(savedir, tag, 'resp_qe_effective.txt')).T

fn = lambda idx: os.path.join(savedir, tag, f'plm_rfid_p_p_{idx:04}.fits')
fn_N0 = lambda idx: os.path.join(savedir, tag, f'Nlzero_semianalytic_rfid_{idx:04}.txt')                  

fn_eff = lambda idx: os.path.join(savedir, tag, f'plm_reff_p_p_{idx:04}.fits')
fn_N0_eff = lambda idx: os.path.join(savedir, tag, f'Nlzero_semianalytic_reff_{idx:04}.txt')   

for idx in np.arange(nsims)[mpi.rank::mpi.size]:
    print(f'rank {mpi.rank} doing {fn(idx)}')
    if not os.path.exists(fn(idx)) or not os.path.exists(fn_eff(idx)):
        plm_i = qlms_dd.get_sim_qlm(k, idx)
        
        plm_rfid = hp.almxfl(plm_i, utils.cli(resp), hp.Alm.getlmax(np.size(plm_i)))
        hp.write_alm(fn(idx), plm_rfid)

        plm_reff = hp.almxfl(plm_i, utils.cli(resp_eff), hp.Alm.getlmax(np.size(plm_i)))
        hp.write_alm(fn_eff(idx), plm_reff)


    if not os.path.exists(fn_N0(idx)) or not os.path.exists(fn_N0_eff(idx)) :
        _n0_semi = param_file.nhl_dd.get_sim_nhl(idx, k, k)
        _N0 = _n0_semi * utils.cli(resp)**2
        _arr = np.array([np.arange(_N0.shape[0]), _N0])
        np.savetxt(fn_N0(idx), _arr.T, header='N0 semi-analytical, normalized with fiducial reponse\n L, N0')

        _N0_eff = _n0_semi * utils.cli(resp_eff)**2
        _arr = np.array([np.arange(_N0_eff.shape[0]), _N0_eff])
        np.savetxt(fn_N0_eff(idx), _arr.T, header='N0 semi-analytical, normalized with effective response\n L, N0')