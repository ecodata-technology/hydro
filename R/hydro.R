
#' Extraterrestrial solar radiation, per Shuttleworth 1993
#' @export
#'
#' @param psi latitude, in degrees
#' @param J Julian day
#' @return solar radiation in mm/day
ETSR <- function(psi, J){

      # Shuttleworth's test:
      # the following should output c(15.0, 15.1, 11.2)
      # round(ETSR(psi=c(30,0,-30), J=105), 1)

      # convert degrees latitude to radians
      psi <- psi * pi / 180

      # solar declination
      delta <- 0.4093 * sin((J*2*pi/365) - 1.405)

      # sunset hour angle
      omega <- acos(-tan(psi) * tan(delta))

      # relative earth-to-sun distance
      dr <- 1 + 0.033 * cos(J*2*pi/365)

      # extraterrestrial solar radiation (mm / day)
      S0 <- 15.392 * dr * (omega * sin(psi) * sin(delta) + cos(psi) * cos(delta) * sin(omega))

      return(S0)
}



#' Mean S0 for a month of the year
#' @export
#'
#' @param month integer from 1 to 12
#' @param latitude decimal degrees
#' @return mean S0 for a month of the year
monthly_S0 <- function(month, latitude){

      # julian days for target month
      j <- rep(1:12, times=c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
      jdays <- which(j==month)

      # average ETSR across all days of the month
      mean(ETSR(latitude, jdays))
}


#' Hargreaves equation for evapotranspiration
#'
#' Per Hargreaves & Samani 1985. (Shuttleworth 1993 seems to incorrectly omit the exponent)
#' @export
#'
#' @param S0 numeric
#' @param tmean degrees C
#' @param tmin degrees C
#' @param tmax degrees C
#' @return evapotranspiration (ETP) in mm/day
hargreaves <- function(S0, tmean, tmin, tmax){
      0.0023 * S0 * (tmax - tmin)^.5 * (tmean + 17.8)
}



#' Annual water balance variables for one site
#' @export
#'
#' @param latitude degrees
#' @param ppt vector of length 12, mm
#' @param tmean vector of length 12, degrees C
#' @param tmax vector of length 12, degrees C
#' @param tmin vector of length 12, degrees C
#' @param annual should annual sums be returned? (default TRUE; if FALSE monthly values for each of the 5 variables will be returned)
#' @return A named vector of hydrologic variables: PPT, PET, AET, CWD, RAR
water_balance <- function(latitude, ppt, tmean, tmax, tmin, annual = TRUE){

      # unclear what this section does and it causes an error with terra
      # if(is.na(tmean[1]) & annual) return(rep(NA, 5))
      # if(is.na(tmean[1]) & !annual) return(rep(NA, 5*12))

      # [note: Wang et al 2012 page 21 use a latitude correction;
      # this note is a placeholder for that calculation if we decide to use it]

      # get S0 values for all 12 months
      S0 <- sapply(1:12, monthly_S0, latitude=latitude)

      # convert S0 from list of rasters to raster stack
      S0 <- rast(S0)

      # monthly water balance variables
      pet <- hargreaves(S0, tmean, tmin, tmax) *
            c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      pet[tmean<0] <- 0 # per Wang et al 2012

      aet <- min(pet, ppt)
      cwd <- pet - ppt
      cwd[cwd<0] <- 0
      rar <- ppt - pet
      rar[rar<0] <- 0

      # annual sums
      if(annual) return(c(PPT=sum(ppt), PET=sum(pet), AET=sum(aet), CWD=sum(cwd), RAR=sum(rar)))
      if(!annual) return(setNames(c(ppt, pet, aet, cwd, rar),
                                  paste0(rep(c("PPE", "PET", "AET", "CWD", "RAR"), each = 12),
                                         rep(1:12, 5))))
}




#' Create a latitude raster
#' @export
#'
#' @param template a raster layer
#' @param already_latlong leave as F unless rasters are already in lat-long projection
#' @return A raster layer with values representing degrees latitude
latitude <- function(template, already_latlong=FALSE){

      requireNamespace("terra")
      requireNamespace("sf")

      template <- template * 1 # force into RAM

      if(!already_latlong) {

            # convert template to new CRS
            template_wgs84 <- terra::project(template, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
            # extract lat values into new raster
            lat <- init(template_wgs84, 'y')
            # set NA values to match template
            lat <- mask(lat, mask=template_wgs84)
            # set CRS back to original (so it matches other input rasters, lat values in the raster are unchanged)
            lat <- terra::project(x=lat, y=template)
      } else {
            lat <- init(template, 'y')
      }

      return(lat)
}


#' Derive annual hydrologic variables
#'
#' This is the main external function in the package, which uses the Hargreaves
#' equation to generate a set of 5 annual water balance variables from a 48
#' monthly temperature and precipitation variables. The function calculates
#' cumulative annual totals for precipitation (PPT), potential
#' evapotranspiration (PET), actual evapotranspiration (AET), climatic water
#' deficit (CWD), and reacharge and runoff (RAR), all in mm. Note that it does
#' not work inside the arctic and antarctic circles.
#' @export
#'
#' @param rasters stack of 48 rasters in this order: ppt1-12, tmean1-12,
#'   tmax1-12, tmin1-12
#' @param temp_scalar multiplier to convert input temperatures to deg C
#' @param ppt_scalar multiplier to convert input precipitation to mm
#' @param ncores number of computing cores to use for parallel processing
#' @param already_latlong leave as F unless rasters are already in lat-long
#'   projection
#' @param annual should annual sums be returned? (default TRUE; if FALSE monthly values for each of the 5 variables will be returned)
#' @param ... additional arguments to clusterR or writeRaster
#' @return A raster stack with 5 layers: PPT, PET, AET, CWD, RAR
#' @references Wang et al 2012 -- ClimateWNA -- Journal of Applied Meteorology
#'   and Climatology. Shuttleworth 1993 -- Evaporation (chapter 4 in Handbook of
#'   Hydrology). Hargreaves and Samani 1985 -- Reference Crop Evaporation from
#'   Ambient Air Temperature
hydro <- function(rasters, temp_scalar=1, ppt_scalar=1, ncores=1, already_latlong=F, annual=T, ...){

      requireNamespace("terra")

      # latitude raster
      lat <- latitude(rasters[[1]], already_latlong=already_latlong)
      rasters <- c(rasters, lat)

      # compute annual water balance variables
      w <- function(x, ...) water_balance(latitude=x[[49]],
                                  ppt=x[[1:12]],
                                  tmean=x[[13:24]],
                                  tmax=x[[25:36]],
                                  tmin=x[[37:48]],
                                  annual = annual)
      if(temp_scalar != 1 | ppt_scalar != 1) w <- function(x, ...) water_balance(latitude=x[[49]],
                                                                         ppt=x[[1:12]] * ppt_scalar,
                                                                         tmean=x[[13:24]] * temp_scalar,
                                                                         tmax=x[[25:36]] * temp_scalar,
                                                                         tmin=x[[37:48]] * temp_scalar,
                                                                         annual = annual)
      # wrap in another terra::rast to return a multi-layer raster instead of list of rasters
      # this matches output format of previous raster implementation
      wb <- rast(w(rasters))

       # remove parallelization (for now). In many cases, single-threaded terra is already faster than parallelized raster.
#
#       beginCluster(ncores, type="SOCK")
#       wb <- clusterR(rasters, calc, args=list(fun=w),
#                      export=c("ETSR", "hargreaves", "water_balance", "monthly_S0"),
#                      ...)
#       endCluster()

      if(annual) names(wb) <- c("PPT", "PET", "AET", "CWD", "RAR")
      if(!annual) names(wb) <- paste0(rep(c("PPT", "PET", "AET", "CWD", "RAR"), each = 12), rep(1:12, 5))
      return(wb)
}


