##
## intsurv: Integrative Survival Models
## Copyright (C) 2017-2019  Wenjie Wang <wjwang.stat@gmail.com>
##
## This file is part of the R package intsurv.
##
## The R package intsurv is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package intsurv is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

## collation after class.R
##' @include class.R
NULL

##' Cox proportional hazard model
##'
##' This function fits regular Cox proportional hazard model (Cox, 1972) by the
##' monotonic quadratic approximation algorithm (Bohning and Lindsay, 1988).  It
##' is not intended to replace \code{\link[survival]{coxph}} or
##' \code{\link[survival]{coxph.fit}} for fitting regular Cox model.  However,
##' it allows non-integer event indicator (between 0 and 1) when computing the
##' partial likelihood (Cox, 1975) and the profiled baseline hazard function,
##' which can be useful for some research problems where event indicators
##' contain uncertainty.
##'
##' @param time A numeric vector for observed times.
##' @param event A numeric vector for event indicators
##' @param x A numeric matrix for design matrix.
##' @param start A numeric vector for starting values.  By default, an all-zeros
##'     vector with appropriate length (equal to the number of columns of the
##'     design matrix) will be used.
##' @param max_iter A positive integer specifying the maximum number of
##'     iterations.
##' @param rel_tol A positive number specifying the convergence tolence on the
##'     relative difference based on L1-norm.
##' @param early_stop A boolean value specifying whether to stop the iteration
##'     once the negative log-likelihood unexpectedly increases, which may
##'     suggest convergence on likelihood, or indicate numerical issues or
##'     implementation bugs.  The default value is \code{TRUE}.
##' @param verbose A boolean value.  If \code{TRUE}, a verbose information will
##'     be given along iterations for tracing the convergence.  The default
##'     value is \code{FALSE}.
##'
##' @references
##'
##' Cox, D. R. (1972). Regression Models and Life-Tables. Journal of the Royal
##' Statistical Society. Series B (Methodological), 34(2), 187--220.
##'
##' Cox, D. R. (1975). Partial Likelihood. Biometrika, 62(2), 269--276.
##'
##' B\"ohning, Dankmar, & Lindsay, B. G. (1988). Monotonicity of Quadratic-
##' Approximation Algorithms. Annals of the Institute of Statistical
##' Mathematics, 40(4), 641--663.
##'
coxph_fit <- function(time, event, x, start = rep(0, ncol(x)), max_iter = 200,
                      rel_tol = 1e-5, early_stop = TRUE, verbose = FALSE)
{
    ## some quick checks on input arguments
    if (! is.matrix(x)) {
        x <- as.matrix(x)
    }
    nObs <- length(time)
    if (length(event) != nObs) {
        stop("The time and event must have the same length.")
    }
    if (nrow(x) != nObs) {
        stop("The length of time must equal ",
             "the number of rows of the design matrix.")
    }
    if (length(start) != ncol(x)) {
        stop("The length of starting values is inappropriate.")
    }
    ## call the routine
    rcpp_coxph(time, event, x, start, max_iter, rel_tol, early_stop, verbose)
}
