##' Wrapping Locations Around the Dateline
##'
##' These functions wrap and unwrap a sequence of longitudes around
##' the dateline.
##'
##' The \code{wrapLon} function wraps the longitudes back into the
##' interval [lmin,lmin+360).  The \code{unwrapLon} function unwraps a
##' sequence of longitudes so the the initial point lies in
##' [lmin,lmin+360), but the subsequent longitudes in the sequence may
##' wander outside that range.
##'
##' @title Dateline adjustment
##' @param lon a vector of longitudes
##' @param lmin western boundary for wrapped longitudes
##' @return a vector of longitudes
##' @export
wrapLon <- function(lon,lmin=-180)
  (lon-lmin)%%360+lmin

##' @rdname wrapLon
##' @export
unwrapLon <- function(lon,lmin=-180)
  cumsum(c(wrapLon(lon[1],lmin),wrapLon(diff(lon))))



##' Correlated Random Walk Filter
##'
##' Fit a correlated random walk to filter a track and predict
##' locations for given time steps.  The movement model is as
##' described in Johnson et al. (2008), but without drift or haul out
##' components, but assumes t distributed errors as described by
##' Albertsen et al. (2015) if `tdf` is positive.
##'
##' The input track must be given as a dataframe where each row is an
##' observed location, with columns
##' \tabular{ll}{
##' \code{date} \tab observation time (as GMT POSIXct) \cr
##' \code{lon} \tab observed longitude \cr
##' \code{lat} \tab observed latitude \cr
##' \code{lon.se} \tab standard error of longitude (optional) \cr
##' \code{lat.se} \tab standard error of latitude (optional) \cr
##' }
##'
##' The filtering model assumes the errors in longitude and latitude
##' are have standard deviations \code{lon.se} and \code{lat.se}
##' scaled by the \eqn{\sigma_{y}}{sigmaY} model parameters. If these
##' columns are missing, they are assumed to be 1.
##'
##' The \code{corPar} and \code{errPar} arguments control how the
##' correlation parameters \eqn{\beta}{beta} and the error scaling
##' parameters \eqn{\sigma_{y}}{sigmaY} apply to the longitudinal
##' and latitudinal processes:
##' \tabular{ll}{
##' \code{"free"} \tab independent parameters are estimated for
##'     longitude and latitude \cr
##' \code{"equal"} \tab a common parameter is estimated for both
##'     longitude and latitude \cr
##' \code{"fixed"} \tab the parameters are determined by \code{par}
##' }
##'
##' @title Correlated Random Walk Filter
##' @param data A dataframe representing the track (see details).
##' @param predict.date Times at which to predict locations (POSIXct).
##' @param par Vector of initial parameter estimates.
##' @param corPar Controls the autocorrelaion parameter for
##'   longitudinal and latitudinal processes.
##' @param errPar Controls the scaling parameter for the observational
##'   errors for the longitudinal and latitudinal processes.
##' @param tdf Degrees of freedom for the multivariate t error
##'   distribution.
##' @param verbose Enable tracing information.
##' @param control List of control parameters for \code{nlminb}.
##' @return Returns a list with components
##'   \item{\code{summary}}{parameter summary table}
##'   \item{\code{par}}{vector of parameter estimates}
##'   \item{\code{track}}{dataframe of the fitted track}
##'   \item{\code{opt}}{the object returned by the optimizer}
##'   \item{\code{tmb}}{the \pkg{TMB} object}
##' The \code{track} dataframe has columns
##'   \item{\code{date}}{time (as GMT POSIXct)}
##'   \item{\code{lon}}{longitude}
##'   \item{\code{lat}}{latitude}
##'   \item{\code{vlon}}{longitudal velocity}
##'   \item{\code{vlat}}{latitudinal velocity}
##'   \item{\code{lon.se}}{standard error of longitude}
##'   \item{\code{lat.se}}{standard error of latitude}
##'   \item{\code{vlon.se}}{standard error of longitudinal velocity}
##'   \item{\code{vlat.se}}{standard error of latitudinal velocity}
##'   \item{\code{observed}}{whether this was an observed time}
##'   \item{\code{predicted}}{whether this was an prediction time}
##' @references
##'     Johnson, D. S., London, J. M., Lea, M. A. and Durban, J. W. (2008).
##'     Continuous-time correlated random walk model for animal telemetry data.
##'     Ecology, 89(5), 1208-1215.
##'
##'     Albertsen, C. M., Whoriskey, K., Yurkowski, D., Nielsen,
##'     A. and Flemming, J. M. (2015).  Fast fitting of non-Gaussian
##'     state-space models to animal movement data via Template Model
##'     Builder.  Ecology, 96(10), 2598-2604.
##'
##'     Lange, K. L., Little, R. J., & Taylor, J. M. (1989). Robust
##'     statistical modeling using the t distribution. Journal of the
##'     American Statistical Association, 84(408), 881-896.
##' @useDynLib RWalc
##' @importFrom TMB MakeADFun sdreport summary.sdreport
##' @importFrom stats nlminb
##' @export
crw <- function(data,predict.date=NULL,
                par=c(1,1,1,1,1,1),
                corPar=c("free","equal","fixed"),
                errPar=c("free","equal","fixed"),
                tdf=-1,verbose=FALSE,
                control = list(eval.max=2000,iter.max=1500,rel.tol=1.0e-3,x.tol=1.5e-2)) {

  corPar <- match.arg(corPar)
  errPar <- match.arg(errPar)

  ## Preprocess data
  data$date <- as.POSIXct(data$date,tz="GMT")
  data$lon <- unwrapLon(data$lon)
  if(is.null(data$fixed)) data$fixed <- FALSE
  if(is.null(data$lon.se)) data$lon.se <- 1
  if(is.null(data$lat.se)) data$lat.se <- 1


  ## Interleave times
  tms <- .POSIXct(sort(union(data$date,predict.date)),tz="GMT")
  obs <- match(data$date,tms)
  prd <- match(predict.date,tms)


  ## Control parameters
  map <- list(
    logBeta=factor(switch(corPar,free=c(1,2),equal=c(1,1),fixed=c(NA,NA))),
    logSigmaY=factor(switch(errPar,free=c(1,2),equal=c(1,1),fixed=c(NA,NA))))


  ## TMB data
  y <- cbind(data$lon,data$lat)
  w <- cbind(data$lon.se,data$lat.se)
  dt <- diff(as.numeric(tms)/60)
  tmb.data <- list(y=y,w=w,dt=dt,obs=obs,tdf=tdf)

  ## TMB parameters
  beta <- par[1:2]
  sigma <- par[3:4]
  sigmaY <- par[5:6]
  mu <- cbind(approx(as.numeric(data$date),data$lon,tms,rule=2)$y,
              approx(as.numeric(data$date),data$lat,tms,rule=2)$y)
  nu <- matrix(0,length(tms),2)
  tmb.pars <- list(logBeta=log(beta),logSigma=log(sigma),logSigmaY=log(sigmaY),mu=mu,nu=nu)

  ## TMB - create objective function
  obj <- MakeADFun(tmb.data,tmb.pars,map,random=c("mu","nu"),DLL="RWalc",silent=!verbose)
  obj$env$inner.control$trace <- verbose
  obj$env$tracemgc <- verbose

  ## Minimize objective function
  opt <- suppressWarnings(
    nlminb(obj$par,obj$fn,obj$gr,control=control))

  ## Extract parameters and track
  sdrep <- sdreport(obj)
  fxd <- summary.sdreport(sdrep,"report")
  rdm <- summary.sdreport(sdrep,"random")
  mu <- matrix(rdm[rownames(rdm)=="mu",1],ncol=2)
  nu <- matrix(rdm[rownames(rdm)=="nu",1],ncol=2)
  mu.se <- matrix(rdm[rownames(rdm)=="mu",2],ncol=2)
  nu.se <- matrix(rdm[rownames(rdm)=="nu",2],ncol=2)

  list(summary=fxd,
       par=fxd[,1],
       track=data.frame(
         date=tms,
         lon=mu[,1],lat=mu[,2],
         vlon=nu[,1],vlat=nu[,2],
         lon.se=mu.se[,1],lat.se=mu.se[,2],
         vlon.se=nu.se[,1],vlat.se=nu.se[,2],
         observed=seq_len(nrow(mu)) %in% obs,
         predicted=seq_len(nrow(mu)) %in% prd),
       opt=opt,obj=obj)
}


##' ARGOS Error Scale Factors for Location Classes
##'
##' Create a dataframe of the multiplicative scaling factors for
##' scaling location accuracy from the reported Argos location class.
##' These are the as used in \pkg{crawl}.
##'
##' @title ARGOS Error Scale Factors
##' @param lc Argos error classes.
##' @return A dataframe with columns
##' \item{\code{lc}}{ARGOS location class}
##' \item{\code{lon.se}}{error scaling factor for longitude}
##' \item{\code{lat.se}}{error scaling factor for latitude}
##' @export
argosScale <- function(lc) {
  lc <- factor(lc,levels=c("3", "2", "1", "0", "A", "B"),ordered=TRUE)
  lon.se <- c(1,1.54,3.72,23.90,13.51,44.22)
  lat.se <- c(1,1.29,2.55,103.70,14.99,32.53)
  data.frame(lon.se=lon.se[lc],lat.se=lat.se[lc])
}



##' State model matrices for the Continuous Time Random Walk Model
##'
##' Constructs the transition matrix \code{A} and innovation
##' covariance matrix \code{Q} for the continuous time random walk
##' model corresponding to parameters \code{beta}, \code{sigma} and
##' time step \code{dt}.
##'
##' @title State model matrices
##' @param beta Parameter vector of length 2.
##' @param sigma Parameter vector of length 2.
##' @param dt Time step.
##' @return A list with components
##'   \item{\code{A}}{transition matrix}
##'   \item{\code{Q}}{innovation covariance matrix}
##' @export
system.matrices <- function(beta,sigma,dt) {
  A <- matrix(0,4,4)
  Q <- matrix(0,4,4)

  A[1,1] <- 1
  A[1,3] <- (1-exp(-beta[1]*dt))/beta[1]
  A[2,2] <- 1
  A[2,4] <- (1-exp(-beta[2]*dt))/beta[2]
  A[3,3] <- exp(-beta[1]*dt)
  A[4,4] <- exp(-beta[2]*dt)

  Q[1,1] <- sigma[1]^2/beta[1]^2*(dt-2*(1-exp(-beta[1]*dt))/beta[1]+(1-exp(-2*beta[1]*dt))/(2*beta[1]))
  Q[2,2] <- sigma[2]^2/beta[2]^2*(dt-2*(1-exp(-beta[2]*dt))/beta[2]+(1-exp(-2*beta[2]*dt))/(2*beta[2]))
  Q[3,3] <- sigma[1]^2/beta[2]^2*beta[2]*(1-exp(-2*beta[1]*dt))/2
  Q[4,4] <- sigma[2]^2/beta[2]^2*beta[2]*(1-exp(-2*beta[2]*dt))/2
  Q[1,3] <- sigma[1]^2/beta[1]^2*(1-2*exp(-beta[1]*dt)+exp(-2*beta[1]*dt))/2
  Q[3,1] <- sigma[1]^2/beta[1]^2*(1-2*exp(-beta[1]*dt)+exp(-2*beta[1]*dt))/2
  Q[2,4] <- sigma[2]^2/beta[2]^2*(1-2*exp(-beta[2]*dt)+exp(-2*beta[2]*dt))/2
  Q[4,2] <- sigma[2]^2/beta[2]^2*(1-2*exp(-beta[2]*dt)+exp(-2*beta[2]*dt))/2

  list(A=A,Q=Q)
}


##' Generate new tracks from a CRW model
##'
##' Given a a template track and the parameters of a crw model this
##' function generates a new track of the same length that coincides
##' with the fitted track at the start point and optionally other
##' specified points along the template track.
##'
##' Locations from the template track can be marked as fixed with the
##' \code{fixed} argument. This should either be \code{NULL} or a
##' vector with an entry for each location, with entries indicating
##' \itemize{
##'    \item{0} the location is not fixed
##'    \item{1} the location is fixed
##'    \item{2} the location is fixed and the point may be a cusp
##'             (autocorrelation is reset).
##' }
##' In the current implementation the first location must always be
##' fixed.  The \code{fixed.err} parameter specifies the covariance of
##' the error in the fixed points, allowing the user to control how
##' acccurately the simulated track reproduces the four components
##' (locations and velocities) of the fixed points.
##'
##' Additional constraints can be placed on the path by rejection
##' sampling through the function \code{point.check}.  This function
##' must accept a time, longitude and latitude and return a boolean
##' indicating whether the point is acceptable.  For example, the
##' track can be constrained to the ocean by supplying a
##' \code{point.check} function that compares the state to a land mask
##' and returns \code{FALSE} for locations on land.
##'
##' Tracks are simulated in the plane.  There is is no polar
##' correction and without a suitable \code{point.check} function,
##' there is nothing prevent the track extending outside the [-90,90]
##' latitudinal limits.
##'
##' @title Regressive bridge sampler
##' @param data A dataframe representing the template track.
##' @param par The model parameters.
##' @param fixed An integer vector indicating which locations in the
##'   template path are to be held fixed.
##' @param fixed.err Covariance matrix for fixed points.
##' @param point.check A function that accepts a time, longitude and
##'   latitude and returns boolean indicating whether the state is
##'   acceptable.
##' @return Returns a dataframe representing the simulated track with
##'   columns
##'   \item{\code{date}}{time as POSIXct}
##'   \item{\code{lon}}{longitude}
##'   \item{\code{lat}}{latitude}
##'   \item{\code{vlon}}{longitudinal velocity}
##'   \item{\code{vlat}}{latitudinal velocity}
##' @importFrom stats rnorm
##' @export
crwSimulate <- function(data,par,fixed=NULL,
                        fixed.err=diag(1.0E-6,4,4),
                        point.check=function(tm,lon,lat) TRUE) {

  beta <- par[1:2]
  sigma <- par[3:4]

  ## First state must be fixed
  fixed <- if(!is.null(fixed)) fixed else rep(FALSE,nrow(data))
  fixed[1] <- TRUE

  ## Calculate system matrices
  ts <- as.POSIXct(data$date,tz="GMT")
  dt <- diff(as.numeric(ts)/60)
  As <- vector("list",length(dt))
  Qs <- vector("list",length(dt))
  for(u in unique(dt)) {
    AQ <- system.matrices(beta,sigma,u)
    k <- which(dt==u)
    As[k] <- AQ[1]
    Qs[k] <- AQ[2]
  }

  ## Times and matrix of states
  xs <- cbind(unname(as.matrix(data[,c("lon","lat")])),0,0)
  n <- nrow(xs)

  ## Prior mean and variance for each state
  ms <- xs
  Vs <- array(fixed.err,c(4,4,n))


  ## Forward pass - generate priors from movement model
  for(k in 2:n)
    if(!fixed[k]) {
      ms[k,] <- As[[k-1]]%*%ms[k-1,]
      Vs[,,k] <- tcrossprod(As[[k-1]]%*%Vs[,,k-1],As[[k-1]])+Qs[[k-1]]
    }

  ## Reverse pass - recursively sample with a Kalman/Regression step
  ## starting from x[k0,]
  sample <- function(k0) {
    x <- ms[k0,] + drop(rnorm(4)%*%chol(Vs[,,k0]))
    xs[k0,] <<- x
    for(k in (k0-1):1) {
      ## Kalman gain
      ## K <- Vs[,,k]%*%t(As[[k]])%*%solve(As[[k]]%*%Vs[,,k]%*%t(As[[k]])+Qs[[k]])
      W <- As[[k]]%*%Vs[,,k]
      K <- crossprod(W,solve(tcrossprod(W,As[[k]])+Qs[[k]]))
      ## Mean, variance update
      mu <- ms[k,] + drop(K%*%(x-As[[k]]%*%ms[k,]))
      ## V <- Vs[,,k] - K%*%As[[k]]%*%Vs[,,k]
      W <- (diag(1,4,4)-K%*%As[[k]])
      V <- tcrossprod(W%*%Vs[,,k],W)+tcrossprod(K%*%Qs[[k]],K)
      R <- chol(V)
      ## point.check/rejection loop
      for(r in 1:100) {
        x <- mu + drop(rnorm(length(mu))%*%R)
        if(fixed[k] || point.check(ts[k],x[1],x[2])) break
        ## If fail, return last fixed point
        if(r==100) return(k0)
      }
      xs[k,] <<- x

      ## Allow discontinuity at a fixed point
      if(fixed[k]==2) {
        x <- ms[k,] + drop(rnorm(4)%*%chol(fixed.err))
        k0 <- k
      }
    }
    ## On success, return 0
    0
  }

  k <- n
  for(i in 1:50) {
    k <- if(i < 25) sample(k) else sample(n)
    if(k==0) return(data.frame(date=ts,lon=xs[,1],lat=xs[,2],vlon=xs[,3],vlat=xs[,4]))
  }
  NULL
}



##' Resample a Track by Linear Interpolation
##'
##' Linearly interpolate in longitude and latitude to resample a track
##' back to a regular time step.
##'
##' The input track is given as a dataframe where each row is an
##' observed location, with (at minimum) columns
##' \tabular{ll}{
##' \code{date} \tab observation time (as GMT POSIXct) \cr
##' \code{lon} \tab observed longitude \cr
##' \code{lat} \tab observed latitude \cr
##' }
##'
##' @title Interpolate Track
##' @param data A dataframe representing a track
##' @param tstep the time step to resample to (in seconds)
##' @return a dataframe with columns
##'   \item{\code{date}}{observation time (as POSIXct)}
##'   \item{\code{lon}}{interpolated longitude}
##'   \item{\code{lat}}{interpolated latitude}
##' @importFrom stats approx
##' @export
interpolateTrack <- function(data,tstep=60*60) {

  if(!is.null(data)) {
    ts <- seq(min(data$date),max(data$date),tstep)
    data.frame(date=ts,
               lon=approx(as.numeric(data$date),data$lon,as.numeric(ts))$y,
               lat=approx(as.numeric(data$date),data$lat,as.numeric(ts))$y)
  }
}







