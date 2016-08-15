# RWalc

[![Travis-CI Build Status](https://travis-ci.org/SWotherspoon/RWalc.svg?branch=master)](https://travis-ci.org/SWotherspoon/RWalc)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/SWotherspoon/RWalc?branch=master&svg=true)](https://ci.appveyor.com/project/SWotherspoon/RWalc)

RWalc provides facilities for fitting an Ornstein-Uhlenbeck process
model to animal tracking data and simulating tracks from the fit.

The package is a bastardization of the
[crawl](https://cran.r-project.org/web/packages/crawl/index.html) and
[argosTrack](https://github.com/calbertsen/argosTrack) packages.  If
you only wish to filter tracks these may be better suited to your
purposes.

## Installing

Installing RWalc requires 

* a recent version of R (> 3.3.0),
* windows users will require [Rtools](https://cran.r-project.org/bin/windows/Rtools/),
* the [TMB](https://cran.r-project.org/web/packages/TMB/index.html)
  package and its dependencies.


RWalc is easily installed from GitHub using the devtools package. 

```R
devtools::install_github("SWotherspoon/RWalc")
```

If you don't have `devtools` installed already, install it first. 

```R
install.packages("devtools")
```

RWalc otherwise does not need devtools for normal use.


## TODO

This is a very early release and any or all features may yet change.
