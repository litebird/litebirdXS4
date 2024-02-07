# noise maps for S4 Chile wide


Simple Healpix nside 2048, in muK, float32 noise maps by jcarron Jan 26 2024.
At the moment 100 per frequencies are available, but more can be generated if needed.

Use at your own risk !

To read them, use

    chwide_maps.noise.get_sim_tmap(freq, idx)
    chwide_maps.noise.get_sim_qumap(freq, idx)

or 

    chwide_maps.noise.get_sim_tqumap(freq, idx)

for I, QU, or TQU.

They contain 
*   white noise matching the depth maps computed by Reijo.
The noise is perfectly independent noise from pixel to pixel, and includes the QU covariance
* Additional non-white noise as a Gaussian field following the PBDR lknee's and exponents, which is then rescaled in real space by the hit map.

