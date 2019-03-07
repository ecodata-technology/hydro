# hydro

This R package computes annual water balance variables based on monthly temperature and precipitation values. The methods are based on those laid out by [Hargreaves and Samani (1985)](https://elibrary.asabe.org/abstract.asp?aid=26773), [Shuttleworth (1993)](https://scholar.google.com/scholar?cluster=5604169752644993837&hl=en&as_sdt=2005&sciodt=0,5), and [Wang et al. (2012)](https://journals.ametsoc.org/doi/full/10.1175/JAMC-D-11-043.1).

The main function `water_balance`, used 48 monthly climate variables (min, mean, and max temperature, and precipitation x 12 months) and latitude to derive 4 integrated hydologic variables: PPT (total precipitation), PET (potential evapotranspiration), AET (actual evapotranspiration), CWD (climatic water deficit), and RNR (recharge and runoff).

The package also includes a function `hydro` that operates on a raster stack of the input variables, returning a stack of the hydrological variables described above.

Currently the package can not calculate values within the Artic and Antartic circles.
