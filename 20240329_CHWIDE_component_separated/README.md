# Documentation for CHWIDE component separated maps

#### Data description
We provide foreground cleaned CMB maps obtained with NILC for CMB-S4 Chile wide survey. There are separate maps for E and B-mode polarizations. 
Currently, about 25 realizations are provided for testing purposes. The masks used for the E and B maps are also provided. All maps are at `NSIDE=2048` 
and `lmax=4096`, with 2.1 arcmin gaussian beam. The current maps are only with the PySM medium complexity foreground model.

#### Input simulations and NILC pipeline
The input maps are simulated with the following components:
* Foregrounds are taken from `combined_foregrounds_mediumcomplexity`. Documentation for these foreground maps can be found [here](https://github.com/CMB-S4/s4mapbasedsims/tree/main/202305_dc0).
* Fetch CMB realizations using: `maps.get_s4_map`.
* Fetch CHWIDE noise maps using: `chwide_map.noise.get_sim_tqumap`.
  
These three components are coadded to produce input simulations for the NILC pipeline (unreleased massively parallel NILC implementation). Some details of the NILC pipeline and limited validation of the NILC E-mode maps can be found [here](https://docs.google.com/presentation/d/1dkzO31pXrOUE63-T2z7nT-9avby-eo7r5J2aDqd642k/edit?usp=drive_link).


#### Data availability
They on the CMB-S4xLiteBIRD shared space on NERSC: `/global/cfs/cdirs/cmbs4xlb/v1/component_separated/chwide`
The E-maps are in `nilc_Emaps/fits` subdir, the B-maps are in `nilc_Bmaps/fits` subdir, the common masks are in `masks_common`.

#### Mask description
There are 3 masks:

1. `chwide_clip0p3relhits_NSIDE2048.fits` &rarr; footprint of the CHWIDE patch after clipping to remove high noise regions.
2. `chwide_clip0p3relhits_3degC2apo_NSIDE2048.fits` &rarr; The above region with 3deg apodization. All E-mode maps use this mask. There is no Galactic mask for the NILC E-maps.
3. `dust-mask-15pc_3pc-apo_NSIDE2048.fits` &rarr; The dust mask masks brightest 15 percent of the sky in polarised dust emission, and an additional 3 percent apodization. This mask is used for the NILC B-maps, in addition to the apodized footprint mask.

For E-maps use mask 2, while for B-maps use mask 2 x mask 3 (ie. an union mask of masks 2 and 3).

#### Map-naming convention
NILC E-mode file naming scheme:  
`NILC_CMB-S4_CHWIDE-E{map_type}_NSIDE2048_fwhm2.1_CHLAT-only_medium_cos-NSIDE2048-lmax4096_mc0{xx}.fits`  
where map_type are strings taking values ‘map’, ‘fg’ and ‘noise’ for NILC cleaned map, residual foreground and noise respectively. The values of xx goes from 00 to 24.  

NILC B-mode file naming scheme:  
`NILC_CMB-S4_CHWIDE-B{map_type}_NSIDE2048_fwhm2.1_CHLAT-only_medium_cos-NSIDE2048-lmax4096_galmask_mc0{xx}.fits`  
where map_type are strings taking values ‘map’, ‘fg’ and ‘noise’ for NILC cleaned map, residual foreground and noise respectively. The values of xx goes from 00 to 24.
