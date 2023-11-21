LiteBIRD map based simulations of noise
=======================================

## Summary
500 realizations of white noise per frequency channel for the LiteBIRD satellite for I, Q, U.
unit = `K_CMB`, nside = 512, coord = Galactic

## Description of simulation
Maps were generated from white noise covariances provided in the first release of the LiteBIRD end-to-end
simulations (`LB_e2e_simulations`). The logic that was followed there was simulating all the detectors
for channels having 48 or less detectors, and simulating roughly 1/3 of the detectors for channels
with more than 48 detectors (focal plane symmetry was preserved).

## Location
Folder location: `/global/cfs/cdirs/cmbs4xlb/v1/noise/lb/{channel}`
where {channel} is one of the 22 channels:
 - LFT: L1-060, L1-078, L2-050, L2-068, L2-089, L3-068, L3-089, L3-119, L4-078, L4-100, L4-140
 - MFT: M1-100, M1-140, M1-195, M2-119, M2-166
 - HFT: H1-195, H1-280, H2-235, H2-337, H3-402

Inside each folder `/global/cfs/cdirs/cmbs4xlb/v1/noise/lb/{channel}` the typical file shape is:
`{channel}_wn_map_0512_mc_{realization number}.fits`
where `{realization number} = 0000, ..., 0499`

How to read the n-th Monte Carlo realization of I, Q, U:
```
sim     = 0        # which mc realization
channel = 'H1-195' # reading channel H1-195
p_lb_noise  = '/global/cfs/cdirs/cmbs4xlb/v1/noise/lb/'
wn_map_path = p_lb_noise+channel+'/'+channel+'_wn_map_0512_mc_'+str(sim).zfill(4)+'.fits'
wn_map, header = hp.read_map(wn_map_path, field=[0,1,2], h=True)
```

There are also covariance matrices of the form `LB_{telescope}_{channel}_covmat_wn.npy`
where `{telescope} = LFT, MFT, HFT`

