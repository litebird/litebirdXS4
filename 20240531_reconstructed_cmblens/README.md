# Reconstucted CMB lensing fields

Louis Legrand (louis.legrand@outlook.com)

The lensing fields are reconstructed from the simulated chwide component separated maps, stored in `/global/cfs/cdirs/cmbs4xlb/v1/component_separated/chwide`, see the [documentation about the component separated Chile Wide map](https://github.com/litebird/litebirdXS4/tree/master/20240329_CHWIDE_component_separated)

The reconstructed lensing plm's are stored in NERSC, in the folder: `/global/cfs/cdirs/cmbs4xlb/v1/lensingrec/`.

The folder `chwide_qe_v1` contains the estimated lensing fields, reconstructed with a quadratic estimator (QE). 
The QE's are computed with the [plancklens code](https://github.com/carronj/plancklens) (using the master branch version as of May 2024) . 


The naming convention is 
`plm_resp_p_p_xxxx.fits` where:
- `plm` for the lensing potential field $\phi_{\rm LM}$
- `resp` means that they have been normalized by the QE response, given in the file `resp_qe.txt` for reference
- `p_p` that it used the (generalized) minimum variance polarization only estimator (combining the E and B maps)
- `xxxx` is the simulation index 

The convergence field can be obtained using $ \kappa_{\rm LM} = -\frac12 L(L+1)\phi_{\rm LM}$. 

The plancklens parameter file `param_chwide_qe.py` contains the implementation of the reconstruction:
- The filtered CMB multipoles used for the lensing reconstruction are from `lmin = 1` to `lmax = 4096`. Be careful that this might create internal delensing bias when delensing the B maps with these lensing maps.
- The same masking is applied to both the E and B maps. It is the union of the apodized chwide footprint `chwide_clip0p3relhits_3degC2apo_NSIDE2048.fits` and the dust mask `dust-mask-15pc_3pc-apo_NSIDE2048.fits`.
- The QE weights are the gradient Cls.
- The filtering noise spectra is assumed to be given by a white noise with $ 1.3 \, \rm \mu K arcmin$ with a beam of $2.1 \, \rm arcmin$. It does not model the 1/f noise.

The plm's are not mean field subtracted. The mean field can be estimated by averaging over different reconstructions. 

For reference, the command used to run the QE is 
``` 
python  $HOME/plancklens/examples/run_qlms.py  $HOME/litebirdXS4/20240531_reconstructed_cmblens/params/param_chwide_qe.py -imin 0 -imax 24 -k p_p -ivp -dd 
```

The jupyter notebook `qe_rec_chwide.ipynb` shows some manipulations of the reconstructed lensing fields, as well as estimates of the noise biases of the QE power spectra (N0 and N1 biases).

There are some plots in the slides `QE_S4xLitebird.pdf`.

## Update v1.1, June 2024

Updated lensing maps are stored in the folder `chwide_qe_v1.1`.
The update concerns the CMB inverse variance filtering scales used in the QE, in particular to avoid internal delensing bias on the B map. The new scale cut is:
- Temperature: lmin=30, lmax=4096
- E modes: lmin=30, lmax=4096
- Bmodes: lmin=200, lmax=4096

The parameter file `param_chwide_qe_lminB200.py` has been used for the plancklens reconstruction. 

## Update August 2024

The folder `chwide_qe_v1.1` now contains the semi analytical N0 bias (normalized) for each simulation. See Planck 2015 gravitational lensing paper, equation A30.
The naming convention is `Nlzero_semianalytic_xxxx.txt`

This folder also contains the fiducial N1 bias, in the file `Nlone_fiducial.txt`.

