from __future__ import annotations
from psutil import cpu_count

import numpy as np
import healpy as hp
from . import unlensed_ffp10

from lenspyx.remapping import deflection, utils_geom
from lenspyx import utils_hp


def build_lensalms(idx, lmax_len, beam_amin, numthreads=0):
    unl_cmbs = unlensed_ffp10.cmb_unl_ffp10(aberration=(0., 0., 0.))
    dlm = unl_cmbs.get_sim_dlm(idx)
    lmax, mmax = unl_cmbs.lmax_unl, unl_cmbs.lmax_unl

    lmaxthingauss = 3 * 5120
    len_geom = utils_geom.Geom.get_thingauss_geometry(lmaxthingauss, 2)
    if numthreads <= 0:
        numthreads = cpu_count(logical=False)
    ffi = deflection(len_geom, dlm, mmax, epsilon=1e-7, numthreads=numthreads)

    mmax_len = lmax_len
    len_tlm =  ffi.lensgclm(unl_cmbs.get_sim_tlm(idx), mmax, 0, lmax_len, mmax_len)
    len_eblm = ffi.lensgclm(unl_cmbs.get_sim_elm(idx), mmax, 2, lmax_len, mmax_len)

    if idx%2 == 0: # Adding tensor modes for odd sims
        len_eblm[0] += utils_hp.alm_copy(unl_cmbs.get_sim_tensor_elm(idx, r=0.003), None, lmax_len, mmax_len)
        len_eblm[1] += utils_hp.alm_copy(unl_cmbs.get_sim_tensor_blm(idx, r=0.003), None, lmax_len, mmax_len)
        len_tlm     += utils_hp.alm_copy(unl_cmbs.get_sim_tensor_tlm(idx, r=0.003), None, lmax_len, mmax_len)

    transf_p = hp.gauss_beam(beam_amin / 180 / 60 * np.pi, lmax=lmax_len, pol=True)
    hp.almxfl(len_tlm, transf_p[:, 0], inplace=True)
    hp.almxfl(len_eblm[0], transf_p[:, 1], inplace=True)
    hp.almxfl(len_eblm[1], transf_p[:, 2], inplace=True)

    return len_tlm, len_eblm


def get_LiteBird_cmbalms(idx, freq:int or str):
    pass
    #return build_lensalms(idx, config.LiteBird_lmaxlen[int(freq)], config.LiteBird_beams[int(freq)])

def get_S4LAT_cmbalms(idx, freq:int or str):
    pass
    #return build_lensalms(idx, config.S4LAT_lmaxlen[int(freq)], config.S4LAT_beams[int(freq)])

def get_S4SAT_cmbalms(idx, freq:int or str):
    pass
    #return build_lensalms(idx, config.S4SAT_lmaxlen[int(freq)], config.S4SAT_beams[int(freq)])
