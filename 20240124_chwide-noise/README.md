# noise maps for S4 Chile wide


Simple Healpix nside 2048 noise maps by jcarron Jan 26 2024. To read them, use

    chwide_maps.noise.get_sim_tmap(freq, idx)
    chwide_maps.noise.get_sim_qumap(freq, idx)

or 

    chwide_maps.noise.get_sim_tqumap(freq, idx)

They contain white noise matching the depth maps computed by Reijo
The noise is perfectly independent noise from pixel to pixel, and includes the QU covariance
Additional noise is included by a Gaussian field following the PBDR lknee's and exponents, which is the rescaled by the hit map.