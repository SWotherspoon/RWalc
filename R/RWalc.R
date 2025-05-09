##' These functions wrap and unwrap a sequence of longitudes around
##' the dateline.
##'
##' The `wrapLon` function wraps the longitudes back into the
##' interval \[lmin,lmin+360).  The `unwrapLon` function unwraps a
##' sequence of longitudes so the initial point lies in
##' \[lmin,lmin+360), but the subsequent longitudes in the sequence may
##' wander outside that range.
##'
##' @title Wrap Locations Around the Dateline.
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


##' Interpolate a track to a given set of time points by one of
##' several methods.
##'
##' The track may consist of several independent segments.  These
##' segments may represent either distinct segments of a single track,
##' or distinct tracks that may overlap in time.
##'
##' The input track must is given as a dataframe where each row is an
##' observed location, with columns:
##'
##' - `segment`: integer label for the segment (optional)
##' - `date`: observation time (GMT POSIXct)
##' - `x`: observed x coordinate
##' - `y`: observed y coordinate
##' - `x.se`: standard error of the x coordinate (optional)
##' - `y.se`: standard error of the y coordinate (optional)
##'
##' It is assumed the input dataframe is ordered by segment and by
##' date within segment.
##'
##' The `predict` argument specifies prediction times for which
##' locations along the track will be predicted.  When the track
##' consists of a single segment, these argument may be a vector of
##' POSIXct times, otherwise it must be a dataframe with columns:
##'
##' - `segment`: track segment (integer, optional)
##' - `date`: prediction time (as GMT POSIXct)
##'
##' The fitted track is returned as a dataframe containing both the
##' original and predicted locations.  To obtain just the predicted
##' locations, the dataframe should be subset by the `predicted`
##' column.
##'
##' Several interpolation/smoothing methods are available:
##'
##' - `"approx"`: linear interpolation in x and y
##' - `"loess"`: loess smoothing in x and y
##' - `"gc"`: Assumes x is longitude and y is latitude and
##'   interpolates along a great circle
##' - `"mean"`: the track is replaced by its weighted centroid
##'
##'
##' @title Track Interpolation
##' @param data A dataframe representing the track (see details).
##' @param predict A vector of times (POSIXct) or a dataframe of
##'   segments and times for which to predict locations.
##' @param method Method used to interpolate the track.
##' @param loess.span Span used in the loess smooth.
##' @return Returns a dataframe with columns:
##'   - `segment`: track segment
##'   - `date`: time (as GMT POSIXct)
##'   - `observed`: whether this was an observed time
##'   - `predicted`: whether this was a predicted time
##'   - `x`: x coordinate
##'   - `y`: y coordinate
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



##' `rwalcControl` selects the numerical minimizer and associated
##' control parameters used by `rwalc`.
##'
##' The numerical minimization function used to fit the model is
##' selected by the `method` argument.  Additional control
##' parameters specific to the chosen minimizer can be set though the
##' dots argument.  See `nlminb` and `optim`
##' for available options.
##'
##' @title Control Values for `rwalc`.
##' @param optim the numerical optimizer used in the fit
##' @param verbose Enable tracing information.
##' @param ... control parameters for the chosen optimizer
##' @return Returns a list with components:
##'   - `optim`: the name of the numerical optimizer as a
##'   string, "nlminb" or "optim"
##'   - `verbose`: should tracing information be reported
##'   - `control`: list of control parameters for the optimizer
##' @seealso `nlminb`, `optim`.
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
##' implemented in the package `crawl`.  Unlike the crawl model,
##' the model implemented here has no drift or haul out components,
##' and assumes t distributed errors as described by Albertsen et
##' al. (2015) if `tdf` is positive.
##'
##' The input track may consist of several independent segments.
##' These may represent either non-overlapping segments of a single
##' track, or distinct tracks that may overlap in time.  The fitted
##' random walk is correlated within a segment, but segments are
##' assumed independent. It is assumed the input dataframe is ordered
##' by segment and by date within segment.
##'
##' The input track to be filtered is supplied as a dataframe
##' (`data`) where each row is an observed location, with columns:
##'
##' - `segment`: track segment (integer, optional)
##' - `date`: observation time (as GMT POSIXct)
##' - `x`: observed x coordinate
##' - `y`: observed y coordinate
##' - `x.se`: standard error of the x coordinate (optional)
##' - `y.se`: standard error of the y coordinate (optional)
##'
##' The filtering model assumes the errors in the spatial coordinates
##' are have standard deviations 'x.se' and 'y.se' scaled by the
##' \eqn{\tau}{tau} model parameters. If these columns are missing,
##' they are assumed to be 1.
##'
##' An estimate of the track is required to initialize the fitting
##' process.  This can be supplied by the user through the
##' `track` argument as a dataframe with the same format returned
##' by `interpolateTrack`.  When this argument is
##' `NULL`, the initial track is generated with
##' `interpolateTrack` from the `data` and `predict`
##' arguments.
##'
##' The `predict` argument specifies prediction times for which
##' locations along the track will be predicted when no initial track
##' is given.  When the track consists of a single segment, these
##' argument may be a vector of POSIXct times, otherwise it must be a
##' dataframe with columns:
##'
##' - `segment`: track segment (integer, optional)
##' - `date`: prediction time (as GMT POSIXct)
##'
##' If an initial track is given, the prediction times are determined
##' from that.
##'
##' The arguments `betaPar`, `sigmaPar` and `tauPar`
##' control how the correlation parameters \eqn{\beta}{beta}, the
##' standard deviations of the innovations for the velocity process
##' \eqn{\sigma}{sigma} and the error scaling parameters
##' \eqn{\tau}{tau} apply to the x and y processes:
##'
##' - `"free"`: independent parameters are estimated for x and y
##' - `"equal"`: a common parameter is estimated for both x and y
##' - `"fixed"`: the parameters are determined by `par`
##'
##' @title Correlated Random Walk Filter
##' @param data A dataframe representing the track (see details).
##' @param predict A vector of times (as POSIXct) or a dataframe of
##'   segments and times for which to predict locations.  Ignored if
##'   `track` is provided.
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
##'   `rwalcControl`)
##' @return Returns a list with components:
##'   - `summary`: parameter summary table
##'   - `par`: vector of parameter estimates
##'   - `track`: dataframe of the fitted track
##'   - `opt`: the object returned by the optimizer
##'   - `tmb`: the `TMB` object
##' The `track` dataframe has columns:
##'   - `segment`: track segment
##'   - `date`: time (as GMT POSIXct)
##'   - `x`: x coordinate
##'   - `y`: y coordinate
##'   - `x.v`: x component of velocity
##'   - `y.v`: y component of velocity
##'   - `x.se`: standard error of x coordinate
##'   - `y.se`: standard error of y coordinate
##'   - `x.v.se`: standard error of x component of velocity
##'   - `y.v.se`: standard error of y component of velocity
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



##' Subset the fitted track to return only those locations that
##' correspond to the prediction times and segments.
##'
##' The current implementation can only extract predicted locations
##' for the times and segments specified in the call to `rwalc`.
##'
##' @title Extract Predicted RWalc Track
##' @param object A fitted object of class "rwalc".
##' @param vel Logical indicating whether estimated velocities should
##'   be returned.
##' @param se Logical indicating whether estimated standard errors
##'   should be returned.
##' @param ... Ignored.
##' @return Returns a dataframe of predicted track locations with columns:
##'   - `segment`: track segment
##'   - `date`: time (as GMT POSIXct)
##'   - `x`: x coordinate
##'   - `y`: y coordinate
##'   - `x.v`: x component of velocity
##'   - `y.v`: y component of velocity
##'   - `x.se`: standard error of x coordinate
##'   - `y.se`: standard error of y coordinate
##'   - `x.v.se`: standard error of x component of velocity
##'   - `y.v.se`: standard error of x component of velocity
##' The velocities are omitted when `vel` is `FALSE` and the
##' standard errors are omitted when `se` is `FALSE`.
##' @export
predict.rwalc <- function(object,vel=FALSE,se=FALSE,...) {
  object$track[object$track$predicted,
               c("segment","date","x","y",
                 if(vel) c("x.v","y.v"),
                 if(se) c("x.se","y.se"),
                 if(vel && se) c("x.v.se","y.v.se"))]
}



##' Subset the fitted track to return only those locations that
##' correspond to times and segments present in the observed track.
##'
##' @title Extract Fitted RWalc Track
##' @param object A fitted object of class "rwalc".
##' @param vel Logical indicating whether estimated velocities
##'   should be returned.
##' @param se Logical indicating whether estimated standard errors
##'   should be returned.
##' @param ... Ignored.
##' @return Returns a dataframe of fitted track locations
##' with columns:
##'   - `segment`: track segment
##'   - `date`: time (as GMT POSIXct)
##'   - `x`: x coordinate
##'   - `y`: y coordinate
##'   - `x.v`: x component of velocity
##'   - `y.v`: y component of velocity
##'   - `x.se`: standard error of x coordinate
##'   - `y.se`: standard error of y coordinate
##'   - `x.v.se`: standard error of x component of velocity
##'   - `y.v.se`: standard error of x component of velocity
##' The velocities are omitted when `vel` is `FALSE` and the
##' standard errors are omitted when `se` is `FALSE`.
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
##' of the plots to display can be selected with the `which`
##' argument.
##'
##' @title Plot a Fitted RWalc Track
##' @param x A fitted object of class "rwalc".
##' @param which Select the plots to display (see details).
##' @param segment Select the segments to display (`NULL` displays all)
##' @param ask if `TRUE`, user is asked before each plot is displayed.
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
##'
##' This function requires that the "x" coordinate of the track
##' represents longitude and the "y" coordinate represents latitude.
##'
##' @title ARGOS Error Scale Factors
##' @param data A dataframe representing the track (see details).
##' @param class Column in data containing the Argos location classes.
##' @return The input dataframe with the appended columns:
##' - `x.se`: error scaling factor for longitude
##' - `y.se`: error scaling factor for latitude
##' @export
argosScale <- function(data,class) {
  lc <- factor(data[,class],levels=c("3", "2", "1", "0", "A", "B"),ordered=TRUE)
  lon.se <- c(1,1.54,3.72,23.90,13.51,44.22)
  lat.se <- c(1,1.29,2.55,103.70,14.99,32.53)
  data$x.se <- lon.se[lc]
  data$y.se <- lat.se[lc]
  data
}




##' Construct the transition matrix `A` and innovation
##' covariance matrix `Q` for the continuous time random walk
##' model corresponding to parameters `beta`, `sigma` and
##' time step `dt`.
##'
##' @title State Model Matrices for an RWalc Model
##' @param beta Parameter vector of length 2.
##' @param sigma Parameter vector of length 2.
##' @param dt Time step.
##' @return A list with components:
##'   - `A`: transition matrix
##'   - `Q`: innovation covariance matrix
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



##' Simulate a new track from the parameters of a fitted `rwalc`
##' model.
##'
##' Given a template track and the parameters of a rwalc model this
##' function generates a new track of the same length that coincides
##' with the fitted track at the start point and optionally other
##' specified points along the template track.
##'
##' Locations from the template track can be marked as fixed with the
##' `fixed` argument. This should either be `NULL` or a
##' vector with an entry for each location, with entries indicating:
##'
##' - 0: the location is not fixed
##' - 1: the location is fixed
##' - 2: the location is fixed and the point may be a cusp
##'             (autocorrelation is reset).
##'
##' In the current implementation the first location in each segment
##' must be fixed.  The `fixed.err` parameter specifies the
##' covariance of the error in the fixed points, allowing the user to
##' control how acccurately the simulated track reproduces the four
##' components (locations and velocities) of the fixed points.
##'
##' Additional constraints can be placed on the path by rejection
##' sampling through the function `point.check`.  This function
##' must accept a time, x and y and return a logical indicating
##' whether the point is acceptable.  For example, the track can be
##' constrained to the ocean by supplying a `point.check`
##' function that compares the state to a land mask and returns
##' `FALSE` for locations on land.
##'
##' The `point.accept` can be used to mark points that should not
##' be checked in the rejection step.  Currently this defaults to the
##' value of `fixed` - so by default fixed points are never
##' checked.
##'
##' Tracks are simulated in the plane.  There is is no polar
##' correction and without a suitable `point.check` function,
##' there is nothing prevent the track extending outside the \[-90,90\]
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
##'   should not be checked with the `point.check` function.
##'   Defaults to the value of `fixed`.
##'
##' @return Returns a dataframe representing the simulated track with
##'   columns:
##'   - `segment`: track segment
##'   - `date`: time (as GMT POSIXct)
##'   - `x`: x coordinate
##'   - `y`: y coordinate
##'   - `x.v`: x component of velocity
##'   - `y.v`: y component of velocity
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
  xs <- if(all(c("x.v","y.v") %in% colnames(data)))
          unname(as.matrix(data[,c("x","y","x.v","y.v")]))
        else
          cbind(unname(as.matrix(data[,c("x","y")])),0,0)
  n <- nrow(xs)

  ## Prior mean and variance for each state
  ms <- xs
  Vs <- fixed.err <- array(fixed.err,c(4,4,n))

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

      if(fixed[k]) {
        k0 <- k
        ## Allow discontinuity at a fixed point
        if(fixed[k]==2)
          x <- ms[k,] + drop(rnorm(4)%*%chol(fixed.err[,,k]))
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





##' Generate a surrogate track from a template track by phase
##' randomization.
##'
##' Given a template track generate a new track of the same length
##' that coincides with the fitted track at the start and end point of
##' each segment.
##'
##' Tracks are generated in the plane, and there are no facilities for
##' placing additional constraints on the track.
##'
##' @title Surrogate tracks by phase randomization.
##' @param data A dataframe representing the template track.
##' @return Returns a dataframe representing the surrogate track with
##'   columns:
##'   - `segment`: track segment
##'   - `date`: time (as GMT POSIXct)
##'   - `x`: x coordinate
##'   - `y`: y coordinate
##' @references
##'   Kantz, H., & Schreiber, T. (2004). Nonlinear time series analysis.
##'   Cambridge university press.
##' @importFrom stats mvfft runif
##' @export
rwalcSurrogate <- function(data) {

  if(is.null(data$segment))
    data$segment <- 1

  surrogate <- function(df) {
    ## Extract track increments
    d <- cbind(diff(as.numeric(df$date)),diff(df$x),diff(df$y))
    ## Creat random phases in connjugate pairs
    p <- double(nrow(d))
    k <- seq_len((length(p)-1)%/%2)
    p[k+1] <- p[length(p)+1-k] <- runif(length(k),-pi,pi)
    ## Random phases and sum increments to rebuild the track
    d <- mvfft(exp(1i*p)*mvfft(d),inverse=TRUE)
    data.frame(
      segment=df$segment,
      date=.POSIXct(cumsum(c(as.numeric(df$date[1]),Re(d[,1]))),"GMT"),
      x=cumsum(c(df$x[1],Re(d[,2]))),
      y=cumsum(c(df$y[1],Re(d[,3]))))
  }

  do.call(rbind,lapply(split(data[,c("segment","data","x","y")],data$segment),surrogate))
}
