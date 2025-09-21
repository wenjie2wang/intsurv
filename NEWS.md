# intsurv 0.3.0

## Major changes

- Simplified the interface of functions `cox_cure()` and `cox_cure_net()`
  through helper functions `cox_cure.control()`, `cox_cure.mstep()`,
  `cox_cure_net.penalty()`.


# intsurv 0.2.2

## New features

- Added AIC to outputs of fitted Cox cure models.
- Added a new argument named `cv_nfolds` to `cox_cure_net()` and
  `cox_cure_net.fit()` for model selection by cross-validation.
- Added new arguments named `surv_offset` and `cure_offset` for optional offset
  terms in survival model component and incidence model component, respectively,
  to `cox_cure()`, `cox_cure.fit()`, `cox_cure_net()`, and `cox_cure_net.fit()`.
  One or more `offset` terms can be also specified via corresponding formula.

## Bug fixes

- Fixed issue caused by zeros in specified `surv_l1_penalty_factor` and
  `cure_l1_penalty_factor` for `cox_cure_net()` and `cox_cure_net.fit()`. Thank
  Sy Han Chiou for reporting this issue.


# intsurv 0.2.1

## Minor changes

- Moved C++ header files to `inst/include/` so that they can be linked by other
  packages.
- Added testing suite with help of the tinytest package.
- Adjusted default maximum number of iterations in M-steps for possibly faster
  convergence in function `cox_cure` and `cox_cure_net`.

## Bug fixes

- Fixed weight indices in function `cIndex`.


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
