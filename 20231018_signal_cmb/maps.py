"""This module builds the CMB signal-only maps from the stored lensed alms and instrument models of s4 and litebird.

    maps in uK

 """

from __future__ import annotations
import os
from os.path import dirname as dirn
import numpy as np
import healpy as hp
from astropy.table import QTable
from astropy import units


# Paths to instrument models:
fn = os.path.abspath(__file__)
path2s4 = os.path.join(
    dirn(dirn(fn)),
    "s4mapbasedsims",
    "202305_dc0",
    "instrument_model",
    "cmbs4_instrument_model.tbl",
)
path2lb = os.path.join(
    dirn(dirn(fn)),
    "litebirdXS4-private",
    "instrument_model_20230614",
    "litebird_instrument_model.tbl",
)
path2cmb = os.path.join(os.environ["CFS"], "cmbs4xlb/v1/cmb", "lcdm_teb_%04d.npy")

try:
    s4 = QTable.read(path2s4, format="ascii.ipac")
    s4.add_index("band")
except:
    print("maps.py: could not read CMB-S4 instrument model at " + path2s4)
    s4 = None
try:
    lb = QTable.read(path2lb, format="ascii.ipac")
    lb.add_index("tag")
except:
    print("maps.py: could not read LiteBird instrument model at " + path2lb)
    lb = None

def _slice_alms(teb, lmax_new):
    """Returns the input teb alms sliced to the new lmax.

    teb(numpy array): input teb alms
    lmax_new(int): new lmax
    """

    lmax = hp.Alm.getlmax(len(teb[0]))
    if lmax_new > lmax:
        raise ValueError('lmax_new must be smaller or equal to lmax')
    elif lmax_new == lmax:
        return teb
    else:
        teb_new = np.zeros((3, hp.Alm.getsize(lmax_new)), dtype=teb.dtype)
        indices_full = hp.Alm.getidx(lmax,*hp.Alm.getlm(lmax_new))
        indices_new = hp.Alm.getidx(lmax_new,*hp.Alm.getlm(lmax_new))
        teb_new[:,indices_new] = teb[:,indices_full]
        return teb_new

def _build_maps(idx, beam_amin: float, nside: int, job="TQU", lmax=4096, coord=None):
    """Returns CMB-S4/LiteBIRD simulated T, Q and U maps.

    idx(int): simulation index
    beam_amin(float): beam width in arcmin
    nside(int): healpix nside of the instrument model
    job(str, optional): one of 'T', 'QU' or 'TQU' (default), for temperature-only, pol-only, or all three maps
    lmax(int, optional): maximum multipole of the alm (default 4096)
    coord(list, optional): list of strings, one of 'C' or 'G' (default), for coordinate system of the output maps

    returns:
        numpy array of shape (ncomp, npix), with ncomp set by 'job' and npix by the nside of the instrument model
        CMB maps in uK

    """
    assert os.path.exists(path2cmb % idx), "cmb alms not found at " + path2cmb % idx
    assert job.upper() in ["T", "QU", "TQU"]
    if coord:
        assert type(coord) == list, "coord should be a list of strings"
        rot = hp.Rotator(coord=coord)

    # Builds pixel and beam:
    bl_T, bl_E, bl_B, _ = hp.gauss_beam(
        beam_amin / 180 / 60 * np.pi, pol=True, lmax=lmax
    ).T

    # Additional window function ?
    # pw_T, pw_P = hp.pixwin(nside, lmax=lmax, pol=True)
    pw_T, pw_P = 1.0, 1.0

    # load lencmbs and apply the transfer function
    teb = _slice_alms(np.load(path2cmb%idx),lmax)
    if coord:
        rot.rotate_alm(teb, inplace=True)
    ncomp = ('T' in job.upper()) + 2 * ('QU' in job.upper())
    maps = np.empty((ncomp, hp.nside2npix(nside)), dtype=float)
    if 'T' in job.upper():
        hp.almxfl(teb[0], bl_T * pw_T, inplace=True)
        maps[0] = hp.alm2map(teb[0], nside)
    if 'QU' in job.upper():
        hp.almxfl(teb[1], bl_E * pw_P, inplace=True)
        hp.almxfl(teb[2], bl_B * pw_P, inplace=True)
        maps['T' in job.upper():] = hp.alm2map_spin(teb[1:], nside, 2, lmax)
    return maps


def get_s4_map(band: str, idx: int, job="TQU"):
    """Returns CMB-S4 simulated T, Q and U maps for the requested channel

    band(str): s4 channel (e.g. 'SAT_f030')
    idx(int): simulation index
    job(str, optional): one of 'T', 'QU' or 'TQU' (default), for temperature-only, pol-only, or all three maps

    returns:
        numpy array of shape (ncomp, npix), with ncomp set by 'job' and npix by the nside of the instrument model
        CMB maps in uK

    """
    assert band in list(s4["band"]), ("possible bands: ", list(s4["band"]))
    beam, nside = s4.loc[band]["fwhm"], s4.loc[band]["nside"]
    assert beam.unit == units.Unit("arcmin")
    return _build_maps(idx, beam.value, nside, job=job)


def get_lb_map(band: str, idx: int, job="TQU", lmax=1024, coord=["C", "G"]):
    """Returns LiteBird simulated T, Q and U maps for the requested channel

    band(str): litebird channel (e.g. 'L1-060')
    idx(int): simulation index
    job(str, optional): one of 'T', 'QU' or 'TQU' (default), for temperature-only, pol-only, or all three maps

    returns:
        numpy array of shape (ncomp, npix), with ncomp set by 'job' and npix by the nside of the instrument model
        CMB maps in uK

    """
    assert band in list(lb["tag"]), ("possible bands: ", list(lb["tag"]))
    beam, nside = lb.loc[band]["fwhm"], lb.loc[band]["nside"]
    assert beam.unit == units.Unit("arcmin")
    return _build_maps(idx, beam.value, nside, job=job, lmax=lmax, coord=coord)
