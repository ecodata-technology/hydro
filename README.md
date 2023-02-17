# hydro

This R package computes water balance variables based on monthly temperature and precipitation. The five computed variables include total precipitation (PPT), potential evapotranspiration (PET), actual evapotranspiration (AET), climatic water deficit (CWD), and recharge and runoff (RAR).

The methods are based on those laid out by [Hargreaves and Samani (1985)](https://elibrary.asabe.org/abstract.asp?aid=26773), [Shuttleworth (1993)](https://scholar.google.com/scholar?cluster=5604169752644993837&hl=en&as_sdt=2005&sciodt=0,5), and [Wang et al. (2012)](https://journals.ametsoc.org/doi/full/10.1175/JAMC-D-11-043.1).

## Installation

You can install the package from GitHub like so:

```
devtools::install_github("matthewkling/hydro")
```

## Usage

The main function, `water_balance`, uses 48 monthly climate variables (precipitation and min, mean, and max temperature, x 12 months) and latitude to derive the five integrated hydologic variables mentioned above. Annual summaries are returned by default, but monthly values can be requested if desired.

The package also includes a function `hydro` that operates on a raster stack of the input variables, returning a stack of the hydrological variables described above. Here's an example:

```
library(hydro)
library(raster)
library(tidyverse)

# download monthly wordlclim data 
monthly <- c("prec", "tmean", "tmax", "tmin") %>%
      map(function(x) getData("worldclim", var = x, res = 10)) %>% 
      stack() %>%
      crop(extent(-125, -70, 20, 60)) 
      
# generate annual water balance variables 
wb <- hydro(monthly, temp_scalar = 0.1, already_latlong = TRUE, annual = TRUE)
```

## Caveats

* Currently the package can not calculate values within the Artic and Antartic circles.
