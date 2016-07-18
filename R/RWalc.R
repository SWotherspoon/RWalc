##' These functions wrap and unwrap a sequence of longitudes around
##' the dateline.
##'
##' The \code{wrapLon} function wraps the longitudes back into the
##' interval [lmin,lmin+360).  The \code{unwrapLon} function unwraps a
##' sequence of longitudes so the the initial point lies in
##' [lmin,lmin+360), but the subsequent longitudes in the sequence may
##' wander outside that range.
##'
##' @title Wrapping Locations Around the Dateline.
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


##' \code{interpolateTrack} interpolates a track to a given set of
##' time points by one of several methods.
##'
##' The track may consist of several independent segments.  These may
##' represent either distinct segments of a single track, or
##' distinct tracks that may overlap in time.
##'
##' The input track must be given as a dataframe where each row is an
##' observed location, with columns
##' \tabular{ll}{
##' \code{segment} \tab integer label for segment (optional) \cr
##' \code{date} \tab observation time (GMT POSIXct) \cr
##' \code{x} \tab observed x coordinate \cr
##' \code{y} \tab observed y coordinate \cr
##' \code{x.se} \tab standard error of the x coordinate (optional) \cr
##' \code{y.se} \tab standard error of the y coordinate (optional) \cr
##' }
##' It is assumed the input dataframe is ordered by segment and by
##' date within segment.
##'
##' The \code{predict} argument specifies prediction times for which
##' locations along the track will be predicted.  When the track
##' consists of a single segment, these argument may be a vector of
##' POSIXct times, otherwise it must be a dataframe with columns
##'
##' \tabular{ll}{
##' segment \tab track segment (integer, optional) \cr
##' date \tab prediction time (as GMT POSIXct)
##' }
##'
##' The fitted track is returned as a dataframe containing both the
##' original and predicted locations.  To obtain just the predicted
##' locations, the dataframe should be subset by the \code{predicted}
##' column.
##'
##' Several interpolation/smoothing methods are available
##' \describe{
##' \item{"approx"}{linear interpolation in x and y}
##' \item{"loess"}{loess smoothing in x and y}
##' \item{"gc"}{Assumes x is longitude and y is latitude and
##' interpolates along a great circle}
##' \item{"mean"}{the track is replaced by its weighted centroid}
##' }
##'
##'
##' @title Track Interpolation
##' @param data A dataframe representing the track (see details).
##' @param predict A vector of times (POSIXct) or a dataframe of
##'   segments and times for which to predict locations.
##' @param method Method used to interpolate the track.
##' @param loess.span Span used in the loess smooth.
##' @return Returns a dataframe with columns
##'   \item{segment}{track segment}
##'   \item{date}{time (as GMT POSIXct)}
##'   \item{observed}{whether this was an observed time}
##'   \item{predicted}{whether this was a predicted time}
##'   \item{x}{x coordinate}
##'   \item{y}{y coordinate}
##' @importFrom stats approx loess predict loess.control weighted.mean
##' @export
interpolateTrack <- function(data,predict=NULL,
                             method=c("approx","loess","gc","mean"),
                             loess.span=0.1) {

  ## Select interpolation method
  method <- match.arg(method)

  ## Preprocess data
  data$date <- as.POSIXct(data$date,tz="GMT")
  if(is.null(data$x.se)) data$x.se <- 1
  if(is.null(data$y.se)) data$y.se <- 1
  if(is.null(data$segment)) data$segment <- 1
  if(is.unsorted(order(data$segment,data$date)))
    warning("Data should be ordered by segment and date within segment")

  ## Convert prediction times to dataframe of dates and segments
  if(is.null(predict)) predict <- data.frame(segment=numeric(0),date=numeric(0))
  if(!is.data.frame(predict))
    predict <- data.frame(
      segment=round(approx(as.numeric(data$date),data$segment,as.numeric(predict),rule=2)$y),
      date=predict)

  ## Interleave times
  track <- unique(rbind(data[,c("segment","date")],predict[,c("segment","date")]))
  track <- track[order(track$segment,track$date),]
  rownames(track) <- NULL

  ## Which locations are observed and which predicted
  tab <- paste(track$segment,as.numeric(track$date),sep="\r")
  track$observed <- tab %in% paste(data$segment,as.numeric(data$date),sep="\r")
  track$predicted <- tab %in% paste(predict$segment,as.numeric(predict$date),sep="\r")

  ## Create interpolation function
  interp <- switch(
    method,
    ## Linear interpolation in x,y
    approx=function(data,tms) {
      cbind(approx(as.numeric(data$date),data$x,as.numeric(tms))$y,
            approx(as.numeric(data$date),data$y,as.numeric(tms))$y)
    },
    ## Loess smooth in x, y
    loess=function(data,tms) {
      fit.x <- loess(x~as.numeric(date),data=data,weights=1/data$x.se^2,
                     span=loess.span,na.action="na.exclude",
                     control=loess.control(surface="direct"))
      fit.y <- loess(y~as.numeric(date),data=data,weights=1/data$y.se^2,
                     span=loess.span,na.action="na.exclude",
                     control=loess.control(surface="direct"))
      cbind(x=predict(fit.x,newdata=data.frame(date=as.numeric(tms))),
            y=predict(fit.y,newdata=data.frame(date=as.numeric(tms))))
    },
    ## Assumes x=lon, y=lat and interpolates along a great circle
    gc=function(data,tms) {
      knots <- as.numeric(data$date)
      ks <- unclass(cut(as.numeric(tms),knots,include.lowest=TRUE))
      f <- (as.numeric(tms)-knots[ks])/(knots[ks+1]-knots[ks])
      p <- pi/180*cbind(data$x,data$y)
      p1 <- p[ks,,drop=FALSE]
      p2 <- p[ks+1,,drop=FALSE]

      d <- acos(sin(p1[,2])*sin(p2[,2])+cos(p1[,2])*cos(p2[,2])*cos(p1[,1]-p2[,1]))
      A <- ifelse(abs(d) < 1.0E-14,1-f,sin((1-f)*d)/sin(d))
      B <- ifelse(abs(d) < 1.0E-14,f,sin(f*d)/sin(d))
      x <- A*cos(p1[,2])*cos(p1[,1])+B*cos(p2[,2])*cos(p2[,1])
      y <- A*cos(p1[,2])*sin(p1[,1])+B*cos(p2[,2])*sin(p2[,1])
      z <- A*sin(p1[,2])+B*sin(p2[,2])
      cbind(x=(180/pi)*atan2(y, x),y=(180/pi)*atan2(z, sqrt(x^2+y^2)))
    },
    ## Replace the track with its centroid
    mean=function(data,tms) {
      cbind(x=rep_len(weighted.mean(data$x,1/data$x.se^2,na.rm=TRUE),length(tms)),
            y=rep_len(weighted.mean(data$y,1/data$y.se^2,na.rm=TRUE),length(tms)))
    })

  ## Interpolate in each segment to generate initial mu
  track$x <- double(nrow(track))
  track$y <- double(nrow(track))
  for(s in unique(track$segment))
    track[track$segment==s,c("x","y")] <- interp(data[data$segment==s,],track$date[track$segment==s])
  track
}




##' \code{rwalcControl} selects the numerical minimizer and associated
##' control parameters used by \code{rwalc}.
##'
##' The numerical minimization function used to fit the model is
##' selected by the \code{method} argument.  Additional control
##' parameters specific to the chosen minimizer can be set though the
##' dots argument.  See \code{\link{nlminb}} and \code{\link{optim}}
##' for available options.
##'
##' @title Control Values for \code{rwalc}.
##' @param optim the numerical optimizer used in the fit
##' @param verbose Enable tracing information.
##' @param ... control parameters for the chosen optimizer
##' @return Returns a list with components
##'   \item{\code{optim}}{the name of the numerical optimizer as a
##'   string, "nlminb" or "optim"}
##'   \item{\code{verbose}}{should tracing information be reported}
##'   \item{\code{control}}{list of control parameters for the optimizer}
##' @seealso \code{\link{nlminb}}, \code{\link{optim}}.
##' @export
rwalcControl <- function(optim=c("nlminb","optim"),verbose=FALSE,...) {
  optim <- match.arg(optim)
  dots <- list(...)
  ## Set default control values
  pars <- switch(optim,
                 nlminb=list(eval.max=3000,iter.max=2000,rel.tol=1.0e-3,x.tol=1.5e-2),
                 optim=list(maxit=2000,reltol=1.0e-3))
  ## Override control parameters
  pars[names(dots)] <- dots
  list(optim=optim,verbose=verbose,control=pars)
}


##' Fit a continuous time correlated random walk to filter a track and
##' predict locations for given time steps.
##'
##' The filter fits a continuous time correlated random walk movement
##' model similar to that described in Johnson et al. (2008) and
##' implemented in the package \pkg{crawl}.  Unlike the crawl model,
##' the model implemented here has no drift or haul out components,
##' and assumes t distributed errors as described by Albertsen et
##' al. (2015) if \code{tdf} is positive.
##'
##' The input track may consist of several independent segments.
##' These may represent either non-overlapping segments of a single
##' track, or distinct tracks that may overlap in time.  The fitted
##' random walk is correlated within a segment, but segments are
##' assumed independent. It is assumed the input dataframe is ordered
##' by segment and by date within segment.
##'
##' The input track to be filtered is supplied as a dataframe
##' (\code{data}) where each row is an observed location, with columns
##' \tabular{ll}{
##' segment \tab track segment (integer, optional) \cr
##' date \tab observation time (as GMT POSIXct) \cr
##' x \tab observed x coordinate \cr
##' y \tab observed y coordinate \cr
##' x.se \tab standard error of the x coordinate (optional) \cr
##' y.se \tab standard error of the y coordinate (optional) \cr
##' }
##' The filtering model assumes the errors in the spatial coordinates
##' are have standard deviations 'x.se' and 'y.se' scaled by the
##' \eqn{\tau}{tau} model parameters. If these columns are missing,
##' they are assumed to be 1.
##'
##' An estimate of the track is required to initialize the fitting
##' process.  This can be supplied by the user through the
##' \code{track} argument as a dataframe with the same format returned
##' by \code{\link{interpolateTrack}}.  When this argument is
##' \code{NULL}, the initial track is generated with
##' \code{interpolateTrack} from the \code{data} and \code{predict}
##' arguments.
##'
##' The \code{predict} argument specifies prediction times for which
##' locations along the track will be predicted when no initial track
##' is given.  When the track consists of a single segment, these
##' argument may be a vector of POSIXct times, otherwise it must be a
##' dataframe with columns
##' \tabular{ll}{
##' segment \tab track segment (integer, optional) \cr
##' date \tab prediction time (as GMT POSIXct)
##' }
##' If an initial track is given, the prediction times are determined
##' from that.
##'
##' The arguments \code{betaPar}, \code{sigmaPar} and \code{tauPar}
##' control how the correlation parameters \eqn{\beta}{beta}, the
##' standard deviations of the innovations for the velocity process
##' \eqn{\sigma}{sigma} and the error scaling parameters
##' \eqn{\tau}{tau} apply to the x and y processes:
##' \tabular{ll}{
##' \code{"free"} \tab independent parameters are estimated for
##'     x and y \cr
##' \code{"equal"} \tab a common parameter is estimated for both
##'     x and y \cr
##' \code{"fixed"} \tab the parameters are determined by \code{par}
##' }
##'
##' @title Correlated Random Walk Filter
##' @param data A dataframe representing the track (see details).
##' @param predict A vector of times (as POSIXct) or a dataframe of
##'   segments and times for which to predict locations.  Ignored if
##'   \code{track} is provided.
##' @param track Dataframe representing an initial estimate of the track (see details).
##' @param par Vector of initial parameter estimates.
##' @param betaPar Controls the autocorrelaion parameter for x and y
##'   processes.
##' @param sigmaPar Controls the standard deviation parameters for the
##'   stochastic innovations of the velocity for the x and y
##'   processes.
##' @param tauPar Controls the scaling parameter for the observational
##'   errors for the x and y processes
##' @param tdf Degrees of freedom for the multivariate t error
##'   distribution.
##' @param bshrink Shrinkage penalty for the correlation parameter.
##' @param control List of control parameters (see
##'   \code{\link{rwalcControl}})
##' @return Returns a list with components
##'   \item{\code{summary}}{parameter summary table}
##'   \item{\code{par}}{vector of parameter estimates}
##'   \item{\code{track}}{dataframe of the fitted track}
##'   \item{\code{opt}}{the object returned by the optimizer}
##'   \item{\code{tmb}}{the \pkg{TMB} object} The \code{track}
##'   dataframe has columns \item{segment}{track segment}
##'   \item{date}{time (as GMT POSIXct)} \item{x}{x coordinate}
##'   \item{y}{y coordinate} \item{x.v}{x component of velocity}
##'   \item{y.v}{y component of velocity} \item{x.se}{standard error
##'   of x coordinate} \item{y.se}{standard error of y coordinate}
##'   \item{x.v.se}{standard error of x component of velocity}
##'   \item{y.v.se}{standard error of y component of velocity}
##'   \item{observed}{whether this was an observed time}
##'   \item{predicted}{whether this was a prediction time}
##' @references
##'   Johnson, D. S., London, J. M., Lea, M. A. and Durban,
##'   J. W. (2008).  Continuous-time correlated random walk model for
##'   animal telemetry data.  Ecology, 89(5), 1208-1215.
##'
##'   Albertsen, C. M., Whoriskey, K., Yurkowski, D., Nielsen, A. and
##'   Flemming, J. M. (2015).  Fast fitting of non-Gaussian
##'   state-space models to animal movement data via Template Model
##'   Builder.  Ecology, 96(10), 2598-2604.
##'
##'   Lange, K. L., Little, R. J., & Taylor, J. M. (1989). Robust
##'   statistical modeling using the t distribution. Journal of the
##'   American Statistical Association, 84(408), 881-896.
##'
##' @useDynLib RWalc
##' @importFrom TMB MakeADFun sdreport summary.sdreport
##' @importFrom stats nlminb optim
##' @export
rwalc <- function(data,
                predict=NULL,
                track=NULL,
                par=c(1,1,1,1,1,1),
                betaPar=c("free","equal","fixed"),
                sigmaPar=c("free","equal","fixed"),
                tauPar=c("free","equal","fixed"),
                tdf=-1,bshrink=1.0E-6,
                control=rwalcControl()) {

  cl <- match.call()

  ## Set parameter constraints
  betaPar <- match.arg(betaPar)
  sigmaPar <- match.arg(sigmaPar)
  tauPar <- match.arg(tauPar)
  map <- list(
    logBeta=factor(switch(betaPar,free=c(1,2),equal=c(1,1),fixed=c(NA,NA))),
    logSigma=factor(switch(sigmaPar,free=c(1,2),equal=c(1,1),fixed=c(NA,NA))),
    logTau=factor(switch(tauPar,free=c(1,2),equal=c(1,1),fixed=c(NA,NA))))

  ## Preprocess data
  data$date <- as.POSIXct(data$date,tz="GMT")
  if(is.null(data$x.se)) data$x.se <- 1
  if(is.null(data$y.se)) data$y.se <- 1
  if(is.null(data$segment)) data$segment <- 1
  if(is.unsorted(order(data$segment,data$date)))
    warning("Data should be ordered by segment and date within segment")

  if(is.null(track)) track <- interpolateTrack(data,predict)
  ## Determine indices of observed locations
  tab <- paste(track$segment,as.numeric(track$date),sep="\r")
  obs <- match(paste(data$segment,as.numeric(data$date),sep="\r"),tab)

  ## TMB data
  y <- cbind(data$x,data$y)
  w <- cbind(data$x.se,data$y.se)
  dt <- diff(as.numeric(track$date)/60)
  seg <- track$segment
  bshrink <- pmax(0,rep_len(bshrink,2))
  tmb.data <- list(y=y,w=w,dt=dt,obs=obs,seg=seg,tdf=tdf,bshrink=bshrink)

  ## TMB parameters
  beta <- par[1:2]
  sigma <- par[3:4]
  tau <- par[5:6]
  mu <- cbind(track$x,track$y)
  nu <- matrix(0,nrow(track),2)
  tmb.pars <- list(logBeta=log(beta),logSigma=log(sigma),logTau=log(tau),mu=mu,nu=nu)

  ## TMB - create objective function
  obj <- MakeADFun(tmb.data,tmb.pars,map,random=c("mu","nu"),DLL="RWalc",silent=!control$verbose)
  obj$env$inner.control$trace <- control$verbose
  obj$env$tracemgc <- control$verbose

  ## Minimize objective function
  opt <- switch(control$optim,
                nlminb=suppressWarnings(nlminb(obj$par,obj$fn,obj$gr,control=control$control)),
                optim=suppressWarnings(do.call(optim,c(obj,control=control["control"]))))

  ## Extract parameters and track
  sdrep <- sdreport(obj)
  fxd <- summary.sdreport(sdrep,"report")
  rdm <- summary.sdreport(sdrep,"random")
  mu <- matrix(rdm[rownames(rdm)=="mu",],ncol=4)
  nu <- matrix(rdm[rownames(rdm)=="nu",],ncol=4)
  track <- cbind.data.frame(
    track[,c("segment","date","observed","predicted")],
    x=mu[,1],y=mu[,2],x.v=nu[,1],y.v=nu[,2],
    x.se=mu[,3],y.se=mu[,4],x.v.se=nu[,3],y.v.se=nu[,4])

  structure(list(call=cl,summary=fxd,par=fxd[,1],track=track,data=data,opt=opt,obj=obj),
            class="rwalc")
}



##' This is a convenience function that subsets the fitted track to
##' return only those locations that correspond to predictions.
##'
##' @title Extract Predicted RWalc Track
##' @param object A fitted object from \code{rwalc}.
##' @param vel Logical indicating whether estimated velocities should
##'   be returned.
##' @param se Logical indicating whether estimated standard errors
##'   should be returned.
##' @param ... Ignored.
##' @return If \code{vel} and \code{se} are true, return a dataframe
##'   of predicted track locations with columns
##'   \item{segment}{track segment}
##'   \item{date}{time (as GMT POSIXct)}
##'   \item{x}{x coordinate}
##'   \item{y}{y coordinate}
##'   \item{x.v}{x component of velocity}
##'   \item{y.v}{y component of velocity}
##'   \item{x.se}{standard error of x coordinate}
##'   \item{y.se}{standard error of y coordinate}
##'   \item{x.v.se}{standard error of x component of velocity}
##'   \item{y.v.se}{standard error of x component of velocity}
##' Otherwise only the appropriate subset of columns are returned
##' @export
##'
predict.rwalc <- function(object,vel=FALSE,se=FALSE,...) {
  object$track[object$track$predicted,
               c("segment","date","x","y",
                 if(vel) c("x.v","y.v"),
                 if(se) c("x.se","y.se"),
                 if(vel && se) c("x.v.se","y.v.se"))]
}



##' This is a convenience function that subsets the fitted track to
##' return only those locations that correspond to the original
##' observations.
##'
##' @title Extract Fitted RWalc Track
##' @param object A fitted object from \code{rwalc}.
##' @param vel Logical indicating whether estimated velocities
##'   should be returned.
##' @param se Logical indicating whether estimated standard errors
##'   should be returned.
##' @param ... Ignored.
##' @return If \code{vel} and \code{se} are true, return a dataframe
##'   of fitted track locations with columns
##'   \item{segment}{track segment}
##'   \item{date}{time (as GMT POSIXct)}
##'   \item{x}{x coordinate}
##'   \item{y}{y coordinate}
##'   \item{x.v}{x component of velocity}
##'   \item{y.v}{y component of velocity}
##'   \item{x.se}{standard error of x coordinate}
##'   \item{y.se}{standard error of y coordinate}
##'   \item{x.v.se}{standard error of x component of velocity}
##'   \item{y.v.se}{standard error of x component of velocity}
##' Otherwise only the appropriate subset of columns are returned.
##' @export
fitted.rwalc <- function(object,vel=FALSE,se=FALSE,...) {
  object$track[object$track$observed,
               c("segment","date","x","y",
                 if(vel) c("x.v","y.v"),
                 if(se) c("x.se","y.se"),
                 if(vel && se) c("x.v.se","y.v.se"))]
}



##' Display the fitted track and observed locations from a fitted RWalc track.
##'
##' Each plot displays the fitted track (blue) and an approximate 95%
##' confidence interval (grey) together with the observed locations
##' (red).  The first two plots display the coordinate profiles of the
##' track over time, while the third plot shows the track.  A subset
##' of the plots to display can be selected with the \code{which}
##' argument.
##'
##' @title Plot a Fitted RWalc Track
##' @param x A fitted object from \code{rwalc}..
##' @param which Select the plots to display (see details).
##' @param segment Select the segments to display (\code{NULL} displays all)
##' @param ask if \code{TRUE}, user is asked before each plot is displayed.
##' @param ... Currently ignored
##' @importFrom grDevices dev.interactive devAskNewPage
##' @importFrom graphics par plot lines points polygon segments
##' @export
plot.rwalc <- function(x,which=1:2,segment=NULL,
                     ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                     ...) {

  plot.profile <- function(date,seg,y,se,lab,date0,y0) {
    se[!is.finite(se)] <- 0
    lwr <- y-2*se
    upr <- y+2*se
    plot(date,y,type="n",ylim=range(c(lwr,upr),na.rm=TRUE),xlab="date",ylab=lab)
    for(s in segment) {
      k <-  which(seg==s)
      polygon(c(date[k],rev(date[k])),c(lwr[k],rev(upr[k])),col="grey90",border=NA)
      lines(date[k],y[k],col="dodgerblue")
    }
    points(date0,y0,col="firebrick",pch=16,cex=0.5)
  }

  plot.track <- function(seg,x,y,x.se,y.se,x0,y0) {
    x.se[!is.finite(x.se)] <- 0
    x.lwr <- x-2*x.se
    x.upr <- x+2*x.se
    y.se[!is.finite(y.se)] <- 0
    y.lwr <- y-2*y.se
    y.upr <- y+2*y.se
    plot(x,y,type="n",xlab="x",ylab="y",
         xlim=range(c(x.lwr,x.upr),na.rm=TRUE),
         ylim=range(c(y.lwr,y.upr),na.rm=TRUE))
    segments(c(x.lwr,x),c(y,y.lwr),c(x.upr,x),c(y,y.upr),col="grey90",lty=1)
    for(s in segment) {
      k <- which(seg==s)
      lines(x[k],y[k],col="dodgerblue")
    }
    points(x0,y0,pch=16,col="firebrick")
  }

  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  if(is.null(segment)) segment <- unique(x$track$segment)
  tr <- x$track[x$track$segment %in% segment,]
  df <- x$data[x$data$segment %in% segment,]
  if(any(which==1))
    plot.profile(tr$date,tr$segment,tr$x,tr$x.se,"x",df$date,df$x)
  if(any(which==2))
    plot.profile(tr$date,tr$segment,tr$y,tr$y.se,"y",df$date,df$y)
  if(any(which==3))
    plot.track(tr$segment,tr$x,tr$y,tr$x.se,tr$y.se,df$x,df$y)
}


##' Create a dataframe of the multiplicative scaling factors for
##' scaling location accuracy from the reported Argos location class.
##' These are the as used in \pkg{crawl}.
##'
##' This function requires that the "x" coordinate represents
##' longitude and the "y" coordinate represents latitude.
##'
##' @title ARGOS Error Scale Factors
##' @param data A dataframe representing the track (see details).
##' @param class Column in data containing the Argos location classes.
##' @return The input dataframe with the appended columns
##' \item{\code{x.se}}{error scaling factor for longitude}
##' \item{\code{y.se}}{error scaling factor for latitude}
##' @export
argosScale <- function(data,class) {
  lc <- factor(data[,class],levels=c("3", "2", "1", "0", "A", "B"),ordered=TRUE)
  lon.se <- c(1,1.54,3.72,23.90,13.51,44.22)
  lat.se <- c(1,1.29,2.55,103.70,14.99,32.53)
  data$x.se <- lon.se[lc]
  data$y.se <- lat.se[lc]
  data
}




##' Constructs the transition matrix \code{A} and innovation
##' covariance matrix \code{Q} for the continuous time random walk
##' model corresponding to parameters \code{beta}, \code{sigma} and
##' time step \code{dt}.
##'
##' @title State Model Matrices for an RWalc Model
##' @param beta Parameter vector of length 2.
##' @param sigma Parameter vector of length 2.
##' @param dt Time step.
##' @return A list with components
##'   \item{\code{A}}{transition matrix}
##'   \item{\code{Q}}{innovation covariance matrix}
##' @export
systemMatrices <- function(beta,sigma,dt) {
  A <- matrix(0,4,4)
  Q <- matrix(0,4,4)
  s <- sigma^2

  A[1,1] <- 1
  A[1,3] <- (1-exp(-beta[1]*dt))/beta[1]
  A[2,2] <- 1
  A[2,4] <- (1-exp(-beta[2]*dt))/beta[2]
  A[3,3] <- exp(-beta[1]*dt)
  A[4,4] <- exp(-beta[2]*dt)

  Q[1,1] <- s[1]*(dt-2*(1-exp(-beta[1]*dt))/beta[1]+(1-exp(-2*beta[1]*dt))/(2*beta[1]))
  Q[2,2] <- s[2]*(dt-2*(1-exp(-beta[2]*dt))/beta[2]+(1-exp(-2*beta[2]*dt))/(2*beta[2]))
  Q[3,3] <- s[1]*beta[1]*(1-exp(-2*beta[1]*dt))/2
  Q[4,4] <- s[2]*beta[2]*(1-exp(-2*beta[2]*dt))/2
  Q[1,3] <- s[1]*(1-2*exp(-beta[1]*dt)+exp(-2*beta[1]*dt))/2
  Q[3,1] <- Q[1,3]
  Q[2,4] <- s[2]*(1-2*exp(-beta[2]*dt)+exp(-2*beta[2]*dt))/2
  Q[4,2] <- Q[2,4]

  list(A=A,Q=Q)
}



##' Given a a template track and the parameters of a rwalc model this
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
##'
##' In the current implementation the first location in each segment
##' must be fixed.  The \code{fixed.err} parameter specifies the
##' covariance of the error in the fixed points, allowing the user to
##' control how acccurately the simulated track reproduces the four
##' components (locations and velocities) of the fixed points.
##'
##' Additional constraints can be placed on the path by rejection
##' sampling through the function \code{point.check}.  This function
##' must accept a time, x and y and return a logical indicating
##' whether the point is acceptable.  For example, the track can be
##' constrained to the ocean by supplying a \code{point.check}
##' function that compares the state to a land mask and returns
##' \code{FALSE} for locations on land.
##'
##' The \code{point.accept} can be used to mark points that should not
##' be checked in the rejection step.  Currently this defaults to the
##' value of \code{fixed} - so by default fixed points are never
##' checked.
##'
##' Tracks are simulated in the plane.  There is is no polar
##' correction and without a suitable \code{point.check} function,
##' there is nothing prevent the track extending outside the [-90,90]
##' latitudinal limits.
##'
##' @title RWalc track sampler
##' @param data A dataframe representing the template track.
##' @param par The model parameters.
##' @param fixed An integer vector indicating which locations in the
##'   template path are to be held fixed.
##' @param fixed.err Covariance matrix for fixed points.
##' @param point.check A function that accepts a time, and an x,y
##'   location and returns a logical indicating whether the location
##'   is acceptable.
##' @param point.accept A logical vector indicating which locations
##'   should not be checked with the \code{point.check} function.
##'   Defaults to the value of \code{fixed}.
##'
##' @return Returns a dataframe representing the simulated track with
##'   columns
##'   \item{segment}{track segment}
##'   \item{date}{time (as GMT POSIXct)}
##'   \item{x}{x coordinate}
##'   \item{y}{y coordinate}
##'   \item{x.v}{x component of velocity}
##'   \item{y.v}{y component of velocity}
##' @importFrom stats rnorm
##' @export
rwalcSimulate <- function(data,par,fixed=NULL,
                        fixed.err=diag(1.0E-6,4,4),
                        point.check=function(tm,x,y) TRUE,
                        point.accept=NULL) {

  beta <- par[1:2]
  sigma <- par[3:4]

  ## First state in each segment must be fixed, and the simulation is
  ## reset at the end of each segment
  seg <- rep_len(if(is.null(data$segment)) 1 else data$segment,nrow(data))
  fixed <- rep_len(if(is.null(fixed)) FALSE else fixed,nrow(data))
  fixed[match(unique(seg),seg)] <- TRUE
  reset <- c(diff(seg)!=0,TRUE)
  ## By default, fixed points are not checked
  point.accept <- rep_len(if(is.null(point.accept)) fixed else point.accept,nrow(data))

  ## Calculate system matrices
  ts <- as.POSIXct(data$date,tz="GMT")
  dt <- diff(as.numeric(ts)/60)
  As <- vector("list",length(dt))
  Qs <- vector("list",length(dt))
  for(u in unique(dt)) {
    AQ <- systemMatrices(beta,sigma,u)
    k <- which(dt==u)
    As[k] <- AQ[1]
    Qs[k] <- AQ[2]
  }

  ## Times and matrix of states
  xs <- if(all(c("x.se","y.se") %in% colnames(data)))
          unname(as.matrix(data[,c("x","y","x.se","y.se")]))
        else
          cbind(unname(as.matrix(data[,c("x","y")])),0,0)
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
    for(k in k0:1) {
      if(reset[k]) {
        k0 <- k
        mu <- ms[k0,]
        R <- chol(Vs[,,k0])
        x <- ms[k0,] + drop(rnorm(4)%*%chol(Vs[,,k0]))
      } else {
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
      }
      ## point.check/rejection loop
      for(r in 1:100) {
        x <- mu + drop(rnorm(length(mu))%*%R)
        if(point.accept[k] || point.check(ts[k],x[1],x[2])) break
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
    if(k==0) {
      df <- cbind.data.frame(segment=seg,date=ts,xs)
      colnames(df) <- c("segment","date","x","y","x.v","y.v")
      return(df)
    }
  }
  NULL
}







