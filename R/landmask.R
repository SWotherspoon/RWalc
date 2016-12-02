##' A land mask based on the GSHHS data set
##'
##' Generate a land mask function based on the Global Self-consistent,
##' Hierarchical, High-resolution Geography Database.  The mask is
##' provided at two (approximate) spatial resolutions: 0.1 degree and
##' 0.05 degrees. The latter requires significantly more memory.
##'
##' Returns a function of three arguments
##' \tabular{lcl}{
##'  \code{tm} \tab \tab time \cr
##'  \code{lon} \tab \tab longitude \cr
##'  \code{lat} \tab \tab latitude \cr
##' }
##' suitable for use with \code{rwalcSimulate}. This mask is actually
##' independent of time and the \code{tm} argument to the mask is
##' ignored.
##'
##' @param res The spatial resolution of the mask, in degrees.
##' @param land The logical value returned for land.
##' @param latmin The minimum allowable latitude.
##' @param latmax The maximum allowable latitude.
##' @return A function of three arguments (tm,lon,lat) that returns
##'   a logical indicating whether (lon,lat) is at sea or on land
##'   and inside the allowable latitude range.
##' @seealso \code{\link{rwalc}}
##' @references
##'   Wessel, P., and W. H. F. Smith, A Global Self-consistent,
##'   Hierarchical, High-resolution Shoreline Database,
##'   J. Geophys. Res., 101, #B4, pp. 8741-8743, 1996.
##'   \url{https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html}
##' @examples
##' mask  <- gshhsMask() ## initialize land mask function
##' mask(0,100,-65) ## test point - time lon,lat
##' @importFrom png readPNG
##' @export
gshhsMask <- function(res=c("0.1","0.05"),land=FALSE,latmin=-90,latmax=90) {
  res <- match.arg(res)

  ## Use only red channel: 0=land, 1=ocean
  mask <- readPNG(system.file("extdata",
                              paste0("land_mask_gshhs-",res,".png"),
                              package="RWalc"))[,,1]
  ## Reflect and convert to binary
  mask <-  if(land) mask[nrow(mask):1,]==0 else mask[nrow(mask):1,]==1

  lonbin <- seq(from=-180,to=180,length.out=ncol(mask)+1)
  latbin <- seq(from=-90,to=90,length.out=nrow(mask)+1)

  function(tm,lon,lat)
    (lat > latmin) &
    (lat < latmax) &
    mask[cbind(.bincode(lat,latbin),.bincode((lon+180)%%360-180,lonbin))]
}


