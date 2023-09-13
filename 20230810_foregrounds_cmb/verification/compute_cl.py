import healpy as hp
from glob import glob
import pickle
import h5py
import toml
import sys

import os.path

from astropy.table import QTable

ID = int(os.environ["SLURM_ARRAY_TASK_ID"])

config = toml.load("common.toml")

chs = QTable.read(
    config["instrument_parameters"],
    format="ascii.ipac",
)

folder = [f for f in glob(config["output_folder"].format(tag="*")) if "C_ell" not in f][ID]
cl = {}
print(folder)
component = os.path.basename(folder)

output_filename = config["output_folder"].format(tag=f"C_ell/C_ell_{component}.pkl")
if os.path.exists(output_filename):
    print(f"{output_filename} already exists")
    sys.exit(0)

for ch in chs:
    try:
        filename = glob(folder + f"/*{ch['tag']}*")[0]
    except IndexError:
        print(folder + f"*{ch['tag']}* NOT FOUND " + ("*" * 20))
        break
    print("reading " + filename)
    try:
        m = hp.read_map(filename, (0, 1, 2))
        nside = hp.npix2nside(len(m[0]))
    except:
        m = hp.read_map(filename)
        nside = hp.npix2nside(len(m))

    # if m.shape[0] == 3 and m[1].sum() == 0:
    #    m = m[0]

    cl[ch["tag"]] = hp.anafast(m, lmax=int(2.5 * nside), use_pixel_weights=False)

if cl:  # empty dicts are false
    with open(output_filename, "wb") as f:
        pickle.dump(cl, f, protocol=-1)
