"""Script exporting lensed CMBs down to lmax 4096 simulations for the litebirdxS4 project


"""
import argparse
import os, glob
import numpy as np
import sys
if os.curdir not in sys.path:
    sys.path.insert(0, os.curdir)
import lencmbs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='cmb export script')
    parser.add_argument('-imin', dest='imin', default=0, type=int, help='starting index)')
    parser.add_argument('-imax', dest='imax', default=-1, type=int, help='last index')
    args = parser.parse_args()
    try:
        from mpi4py import MPI
        rank = MPI.COMM_WORLD.Get_rank()
        size = MPI.COMM_WORLD.Get_size()
        barrier = MPI.COMM_WORLD.Barrier
        finalize = MPI.Finalize
        print('mpi.py : rank %s in %s' % (rank, size))
    except:
        rank, size,  barrier, finalize = 0, 1, lambda: -1, lambda: -1
        print('mpi.py: unable to import mpi4py\n')
    for idx in range(args.imin, args.imax + 1)[rank::size]:
        fn = os.path.join('global/cfs/cdirs/cmbs4xlb/v1/cmb', 'lcdm_teb_%03d.npy'%idx)
        if not os.path.exists(fn) and (0 <= idx <= 499):
            t, eb = lencmbs.build_lensalms(idx, 4096, 0.)
            np.save(np.array([t, eb[0], eb[1]]))
    barrier()
    if rank == 0:
        fns = glob.glob(os.path.join('global/cfs/cdirs/cmbs4xlb/v1/cmb', 'lcdm_teb_???.npy'))
        print('There are %s arrays on disk'%len(fns))
    finalize()