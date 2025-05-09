% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RWalc.R
\name{fitted.rwalc}
\alias{fitted.rwalc}
\title{Extract Fitted RWalc Track}
\usage{
\method{fitted}{rwalc}(object, vel = FALSE, se = FALSE, ...)
}
\arguments{
\item{object}{A fitted object of class "rwalc".}

\item{vel}{Logical indicating whether estimated velocities
should be returned.}

\item{se}{Logical indicating whether estimated standard errors
should be returned.}

\item{...}{Ignored.}
}
\value{
Returns a dataframe of fitted track locations
with columns:
  - `segment`: track segment
  - `date`: time (as GMT POSIXct)
  - `x`: x coordinate
  - `y`: y coordinate
  - `x.v`: x component of velocity
  - `y.v`: y component of velocity
  - `x.se`: standard error of x coordinate
  - `y.se`: standard error of y coordinate
  - `x.v.se`: standard error of x component of velocity
  - `y.v.se`: standard error of x component of velocity
The velocities are omitted when `vel` is `FALSE` and the
standard errors are omitted when `se` is `FALSE`.
}
\description{
Subset the fitted track to return only those locations that
correspond to times and segments present in the observed track.
}
