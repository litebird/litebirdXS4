LiteBIRD external mass tracers
=======================================

## Summary
We provide 500 realizations of the harmonic coefficients of the external mass tracers (signal/noise). 
For a reference, we also save their signal auto/cross spectra which are used to generate the harmonic coefficients. 

## Description of simulation
For the LiteBIRD delensing study, in addition to the S4 CMB lensing map, we plan to explore the importance of the external mass tracers, i.e., 
the CIB and galaxies. In this project, we follow the strategy of https://arxiv.org/abs/2110.09730 where we simply generate these external tracers 
as Gaussian random fields drawn from auto and cross-angular power spectra between CMB-lensing, CIB, and galaxies, 
with the Gaussian kappa in /global/cfs/cdirs/cmbs4xlb/v1/cmb/. 

 * We assume a single CIB map that mimics the Planck GNILC CIB map over ~ 40% of the sky.
   Thus, our simulated CIB is simply a single map. The model is taken from Appendix D of https://arxiv.org/abs/1502.01591 (or https://arxiv.org/abs/1705.02332).
   The residual foreground in the CIB map is added as a random Gaussian field following http://arxiv.org/abs/1705.02332.
   This CIB tracer is used only for the LiteBIRD delensing purpose.
   Note that our simulated CIB map differs from the PySM CIB map which is based on Websky, e.g., in terms of the lack of non-Gaussianity. 
 * The full-sky number density fluctuations of galaxies per redshift bin are prepared assuming LSST-like and Euclid-like galaxy surveys.
   Specifications for the redshift distribution of galaxies, number density of galaxies per arcminutes squared,
   and galaxy linear bias are taken from http://arxiv.org/abs/1705.02332 for LSST, and https://arxiv.org/abs/1606.00180 for Euclid.
   The number of redshift bins for each surveys is chosen as five, where the number density of galaxies are equal for all bins,
   which is optimal for the LiteBIRD delensing found by an investigation of the LiteBIRD delensing forecast.
   The shot noise as a random Gaussian field is added to the map at each redshift bin. 


## Location of the harmonic coefficients
The harmonic coefficients are stored at `/global/cfs/cdirs/cmbs4xlb/v1/mass/alm/`. 
The filename is determined as `{signal or noise}_{mass tracer name}_{realization}.pkl` where 
`{signal or noise}` is signal (s) or noise (n) components of the mass tracer, `{mass tracer name}` specifies the mass tracer: 

 * `ks4`: S4 CMB lensing map
 * `klb`: LiteBIRD CMB lensing map
 * `cib`: CIB
 * `euc{zi}n5`: Euclid-like galaxy density fluctuations at the {zi}-th bin of the photo-z
 * `lss{zi}n5`: LSST-like galaxy density fluctuations at the {zi}-th bin of the photo-z

The `{realization number}` is 0000, ..., 0499. 
Currently, we generate the S4 lensing map as a sum of the input lensing convergence and gaussian random noise whose spectra mimic the S4 reconstruction noise. 
This lensing map is provided just for the debugging purpose. 
LiteBIRD CMB lensing map generated in the same way as the S4 lensing map is also provided for a reference. 

To read these data, please import pickle and use e.g. `pickle.load(open({path to file},"rb"))`. 
Note that this file format is requested by the current delensing pipeline of LiteBIRD, which partially uses Fortran code in addition to Python. 
The pipeline assumes a file in the pickle format that contains the array of the harmonic coefficients in the Healpix format rather than the Healpy one. 
In this format, the harmonic coefficients are given as a 2D array where the first and second indexes are ell (0<=ell<=lmax) and m (0<=m<=lmax), respectively. 
For a given ell and m, the corresponding index for array in the healpy is obtained with healpy.sphtfunc.Alm.getidx(lmax,l,m). 

## Location of the spectra

The input angular auto/cross spectra are saved at `/global/cfs/cdirs/cmbs4xlb/v1/mass/spec/`. 
The filename is determined as cl{mass tracer name 1}{mass tracer name 2}.dat where the first and second column contains multipole and angular power spectrum, respectively. 

To read these files, please use e.g., `l, cl = np.loadtxt({path to file})`. 




