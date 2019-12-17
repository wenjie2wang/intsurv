# intsurv 0.2.1

## Minor changes

- Moved C++ header files to `inst/include/` so that they can be linked by other
  packages.
- Added testing suite with help of the tinytest package.
- Adjusted default maximum number of iterations in M-steps for possibly faster
  convergence in `cox_cure()` and `cox_cure_net()`.

## Bug fixes

- Fixed weight indices in function `cIndex()`.


# intsurv 0.2.0

## New features

- Added function `cox_cure` for fitting regular Cox cure rate model and an
  extended Cox cure rate model for right-censored data with uncertain event
  status.
- Added function `cox_cure_net` for fitting regularized Cox cure rate model for
  right-censored data with possible uncertain event status.
- Added function `cIndex` for computing concordance index with possible weights.
- Added function `simData4cure` for simulating survival data with possible
  uncertain event status.


# intsurv 0.1.0-alpha.2

## Major changes

- Added methods for parameters' initialization.


# intsurv 0.1.0-alpha.1

## New features

- Alpha version of the package for paper submission and ease of reproducing main
  simulation results.
