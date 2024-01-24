import healpy as hp
import numpy as np
import os
from lenspyx.utils_hp import Alm, almxfl
from lenspyx.remapping.utils_geom import Geom
from scipy.interpolate import UnivariateSpline as spl
from multiprocessing import cpu_count
from numpy.random import default_rng

# --- freqs on disk for alt1 chile wide ---
freqs = [30, 40, 90, 150, 220, 280]
XYs = ['II', 'IQ', 'IU', 'QQ', 'QU', 'UU']
path2cov = '/global/cfs/cdirs/cmbs4/AoA/August2022/scaled_outputs/alternative_1/chlat_cd_wide/cov_%03d.fits' # in K^2, nside 128
path2hits = '/global/cfs/cdirs/cmbs4/AoA/August2022/scaled_outputs/alternative_1/chlat_cd_wide/hits_%03d.fits'

nside = 128

# --- PBDR 1/f noise parameters ---
freqs_PBDR = [25, 40, 90, 150, 230, 280]
lknee_T_PBDR = [415, 391, 1932, 3017, 6740, 6792]
lknee_P_PBDR = [700] * len(freqs_PBDR)
a_T = [3.5] * len(freqs_PBDR)
a_P = [1.4] * len(freqs_PBDR)
nlev_T_PBDR = [27.1, 11.6, 2.0, 2.0, 6.9, 16.9]
nlev_P_PBDR = [37.6, 15.5, 2.7, 3.0, 9.8, 23.9]
beam_PBDR = [7.8, 5.3, 2.2, 1.4, 1.0, 0.9]

def get_lknee_a(TorP, freq):
    lknee, a = (lknee_T_PBDR, a_T) if TorP == 'T' else (lknee_P_PBDR, a_P)
    return spl(freqs_PBDR, lknee, k=1, s=0)(freq),spl(freqs_PBDR, a, k=1, s=0)(freq)

def get_nlev_ukamin(freq: int, xx: str, thres=10, nside_out=128):
    """Get's noise level from Cl_noise = mean (pix var map)


    """
    assert xx.upper() in ['II', 'QQ', 'UU']
    assert freq in freqs, (freq, freqs)
    path = path2cov%freq
    cov =  hp.read_map(path, field=XYs.index(xx.upper()))
    assert cov.size == 12 * nside ** 2, (cov.size, 12 * nside ** 2)
    if nside_out != nside:
        cov = hp.ud_grade(cov, nside_out) * (nside_out / nside) ** 2
    cov = cov[np.where(cov > 0)]
    nlev = np.sqrt(np.mean(cov[np.where(cov < thres ** 2 * np.min(cov))]))
    nlev_uKamin = nlev * 1e6 * np.sqrt(hp.nside2pixarea(nside_out, degrees=True) * 60 * 60)
    return nlev_uKamin

def get_rhits(freq: int, nside_out: int):
    assert freq in freqs
    hits =  hp.read_map(path2hits%freq)
    return hp.ud_grade(hits / np.max(hits), nside_out)

def get_cov(freq: int, xx:str, nside_out=2048):
    assert xx.upper() in ['II', 'QQ', 'UU']
    assert freq in freqs, (freq, freqs)
    path = path2cov%freq
    cov =  hp.read_map(path, field=XYs.index(xx.upper())) * 1e12 # in uK^2
    assert cov.size == 12 * nside ** 2, (cov.size, 12 * nside ** 2)
    return hp.ud_grade(cov, nside_out) * (nside_out / nside) ** 2

def _PBDRlike_noise_cls(freq, lmax=5120, _beamdeconvolved=True):
    lknee_I, alpha_I = get_lknee_a('T', freq)
    lknee_P, alpha_P = get_lknee_a('P', freq)
    cl_w = np.ones(lmax + 1, dtype=float)
    cl_knee_I = np.ones(lmax + 1, dtype=float)  * (np.arange(1, lmax + 2) / lknee_I) ** (-alpha_I)
    cl_knee_P = np.ones(lmax + 1, dtype=float)  * (np.arange(1, lmax + 2) / lknee_P) ** (-alpha_P)
    idf = freqs.index(freq)
    bl2 = hp.gauss_beam(beam_PBDR[idf] / 180 / 60 * np.pi, lmax=lmax) ** 2 if _beamdeconvolved else 1.
    return  ((nlev_T_PBDR[idf] / 180 / 60 * np.pi) ** 2 * (cl_w + cl_knee_I) / bl2,
             (nlev_P_PBDR[idf] / 180 / 60 * np.pi) ** 2 * (cl_w + cl_knee_P) / bl2 * (np.arange(lmax + 1) > 1))

def _syn_alm(rng, cl:np.ndarray):
    lmax = cl.size - 1
    mmax = lmax
    rlm_dtype = np.float64

    alm_size = Alm.getsize(lmax, mmax)
    alm = 1j * rng.standard_normal(alm_size, dtype=rlm_dtype)
    alm += rng.standard_normal(alm_size, dtype=rlm_dtype)
    almxfl(alm, np.sqrt(cl * 0.5), mmax, True)
    real_idcs = Alm.getidx(lmax, np.arange(lmax + 1, dtype=int), 0)
    alm[real_idcs] = alm[real_idcs].real * np.sqrt(2.)
    return alm

def get_sim_tmap(freq, lmax=5120, nside_out=2048, seed=None, _facknee=1., _facwhite=1.):

    # ---- sim params and geom ----
    mmax = lmax
    lknee_I, alpha_I = get_lknee_a('T', freq)
    lknee_P, alpha_P = get_lknee_a('P', freq)

    cl_white = np.ones(lmax + 1, dtype=float) * _facwhite
    s2_white_I = np.sum(cl_white * (2 * np.arange(lmax + 1) + 1)) / (4 * np.pi)
    # FIXME: check this and factor root 2
    s2_white_P = np.sum(cl_white[2:] * (2 * np.arange(2, lmax + 1) + 1)) / (4 * np.pi)

    rescal_I = np.sqrt(get_cov(freq, 'II', nside_out=nside_out) / s2_white_I) # Such that white part matches Reijo depth
    rescal_Q = np.sqrt(2.) * np.sqrt(get_cov(freq, 'QQ', nside_out=nside_out) / s2_white_P)
    rescal_U = np.sqrt(2.) * np.sqrt(get_cov(freq, 'UU', nside_out=nside_out) / s2_white_P)

    geom = Geom.get_healpix_geometry(nside_out)
    nthreads = int(os.environ.get('OMP_NUM_THREADS', cpu_count()))

    # ---- syn alm ----
    rng = default_rng(seed=seed)
    cl_knee_I = np.ones(lmax + 1, dtype=float)  * (np.arange(1, lmax + 2) / lknee_I) ** (-alpha_I) * _facknee
    cl_knee_P = np.ones(lmax + 1, dtype=float)  * (np.arange(1, lmax + 2) / lknee_P) ** (-alpha_P) * _facknee

    tlm = _syn_alm(rng, cl_white + cl_knee_I)
    eblm = np.zeros((2, Alm.getsize(lmax, lmax)), dtype=complex)
    eblm[0] = _syn_alm(rng, (cl_white + cl_knee_P) * (np.arange(lmax + 1) >= 2))
    eblm[1] = _syn_alm(rng, (cl_white + cl_knee_P) * (np.arange(lmax + 1) >= 2))

    # ---- maps ----
    tmap = geom.alm2map(tlm, lmax=lmax, mmax=mmax, nthreads=nthreads)
    tmap *= rescal_I

    qmap, umap = geom.alm2map_spin(eblm, 2, lmax, mmax, nthreads=nthreads)
    qmap *= rescal_Q
    umap *= rescal_U
    return tmap, qmap, umap

