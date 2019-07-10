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

##' Cox Proportional Hazard Model
##'
##' Fit a regular Cox proportional hazard model (Cox, 1972) for right-censored
##' data by the monotonic quadratic approximation algorithm (Bohning and
##' Lindsay, 1988).  This function is not intended to replace \code{coxph} or
##' \code{coxph.fit} in \pkg{survival} package for fitting Cox model.  However,
##' it allows non-integer event indicator (between 0 and 1) when computing the
##' partial likelihood (Cox, 1975) and the profiled baseline hazard function,
##' which can be useful for some research problems where event indicators
##' contain uncertainty.  The Breslow's method is used for tied survival times.
##'
##' @usage
##' coxph_fit(time, event, x, offset = NULL, start = rep(0, ncol(x)),
##'           max_iter = 200, rel_tol = 1e-5, early_stop = TRUE,
##'           verbose = FALSE, ...)
##'
##' @param time A numeric vector for observed times.
##' @param event A numeric vector for event indicators
##' @param x A numeric matrix for design matrix.
##' @param offset An opotional offset term added to the linear predictor with
##'     fixed coefficient one.  The default value is \code{NULL} for no offset.
##'     The length of the specified offset must be either one or equal to the
##'     number of observations.
##' @param start An optional numeric vector for starting values.  By default, an
##'     all-zeros vector with appropriate length (equal to the number of columns
##'     of the design matrix) will be used.
##' @param max_iter A positive integer specifying the maximum number of
##'     iterations.
##' @param rel_tol A positive number specifying the tolence that determines the
##'     convergence of the covariate coefficient estimates.  The tolerance is
##'     compared with the relative change between estimates from two consecutive
##'     iterations that is measured by ratio of the L1-norm of their difference
##'     to the sum of their L1-norm.  The default value is \code{1e-5}.
##' @param early_stop A boolean value specifying whether to stop the iteration
##'     once the negative log-likelihood unexpectedly increases, which may
##'     suggest convergence on likelihood, or indicate numerical issues or
##'     implementation bugs.  The default value is \code{TRUE}.
##' @param verbose A boolean value.  If \code{TRUE}, a verbose information will
##'     be given along iterations for tracing the convergence.  The default
##'     value is \code{FALSE}.
##' @param ... Other arguments for future usage.
##'
##' @return
##'
##'
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
##' @export
coxph_fit <- function(time, event, x, offset = NULL,
                      start = rep(0, ncol(x)),
                      max_iter = 200, rel_tol = 1e-5,
                      early_stop = TRUE, verbose = FALSE,
                      ...)
{
    ## warning on `...`
    warn_dots(..., .fun_name = "coxph_fit")

    ## some quick checks on input arguments
    if (! is.matrix(x)) {
        x <- as.matrix(x)
    }
    nObs <- length(time)
    if (length(event) != nObs) {
        stop("The time and event must have the same length.", call. = FALSE)
    }
    if (nrow(x) != nObs) {
        stop("The length of time must equal ",
             "the number of rows of the design matrix.", call. = FALSE)
    }
    if (length(start) != ncol(x)) {
        stop("The length of starting values is inappropriate.", call. = FALSE)
    }
    ## prepare offset
    if (is.null(offset)) {
        offset <- 0
    } else if (length(offset) == 1L && offset != 0) {
        offset <- rep(offset, nObs)
    } else if (length(offset) != nObs) {
        stop("The length of offset term must be ",
             "either one or equal to number of observations.",
             call. = FALSE)
    }
    ## call the routine
    out <- rcpp_coxph(time, event, x, offset, start,
                      max_iter, rel_tol, early_stop, verbose)
    ## return
    class(out) <- "coxph_fit"
    out
}


##' Regularized Cox Proportional Hazard Model
##'
##' Fit regularized Cox model via penalized maximum likelihood.
##'
##'
##'
coxnet_fit <- function(time, event, x, offset = NULL,
                       l1_lambda = 0, l2_lambda = 0,
                       l1_penalty_factor = NULL,
                       start = rep(0, ncol(x)),
                       max_iter = 200, rel_tol = 1e-5,
                       early_stop = TRUE, verbose = FALSE,
                       ...)
{
    ## warning on `...`
    warn_dots(..., .fun_name = "coxnet_fit")



}
