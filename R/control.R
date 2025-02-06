##
## intsurv: Integrative Survival Models
## Copyright (C) 2017-2025  Wenjie Wang <wang@wwenjie.org>
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

##' Common Control Parameters
##'
##' Provides an interface to specify some common control parameters for
##' functions in intsurv.
##'
##' @param maxit A positive integer specifying the maximum iteration number.
##'     The default value is \code{1000}.
##' @param epsilon A positive number specifying the tolerance that determines
##'     the convergence of the coefficient estimates.  The tolerance is compared
##'     with the relative change between estimates from two consecutive
##'     iterations, which is measured by the ratio of the L1-norm of their
##'     difference to the sum of their L1-norms plus one.
##' @param standardize A logical value specifying if each covariate should be
##'     standardized to have mean zero and standard deviation one internally for
##'     numerically stability and fair regularization.  The default value is
##'     \code{TRUE}.  The coefficient estimates will always be returned in
##'     original scales.
##' @param start A numeric vector representing the initial values for the
##'     underlying model estimation procedure.  If \code{standardize} is
##'     \code{TRUE}, the specified initial values will be scaled internally to
##'     match the standardized data.  The default initial values depend on the
##'     specific models and are usually zeros.  If inappropriate initial values
##'     (in terms of length) are specified, the default values will be used.
##' @param offset A numeric vector specifying the offset term.  The length of
##'     the specified offset term should be equal to the sample size.
##' @param verbose A nonnegative integer for verbose outputs, which is mainly
##'     useful for debugging.
##' @param save_call A logical value indicating if the function call should be
##'     saved.  For large datasets, saving the function call would increase the
##'     size of the returned object dramatically.  We may want to set
##'     \code{save_call = FALSE} if the original function call is not needed.
##' @param ... Other arguments.  For \code{cox_cure.control}, they will be
##'     passed to \code{intsurv.control}; for \code{cox_cure_net.control}, they
##'     will be passed to \code{intsurv.control} and \code{cox_cure.control}.
##'     For \code{intsurv.control}, they will be ignored.
##'
##' @return An S3 object of class \code{"intsurv.control"},
##'     \code{cox_cure.control}, or \code{cox_cure_net.control}.
##'
##' @export
intsurv.control <- function(maxit = 300,
                            epsilon = 1e-4,
                            standardize = TRUE,
                            start = NULL,
                            offset = NULL,
                            verbose = 0,
                            save_call = TRUE,
                            ...)
{
    ## TODO add checks
    structure(list(maxit = maxit,
                   epsilon = epsilon,
                   standardize = standardize,
                   start = null2num0(start),
                   offset = null2num0(offset),
                   verbose = verbose,
                   save_call = save_call),
              class = "intsurv.control")
}

##' @rdname intsurv.control
##'
##' @param tail_completion A character string specifying the tail completion
##'     method for conditional survival function.  The available methods are
##'     \code{"zero"} for zero-tail completion after the largest event times (Sy
##'     and Taylor, 2000), \code{"exp"} for exponential-tail completion (Peng,
##'     2003), and \code{"tau-zero"} for zero-tail completion after a specified
##'     \code{tail_tau}.  The default method is the zero-tail completion
##'     proposed by Sy and Taylor (2000).
##' @param tail_tau A numeric number specifying the time of zero-tail
##'     completion.  It will be used only if \code{tail_completion =
##'     "zero-tau"}.  A reasonable choice must be a time point between the
##'     largest event time and the largest survival time.
##' @param pmin A positive number specifying the minimum value of probabilities
##'     for numerical stability.  The default value is \code{1e-5}.
##'
##' @export
cox_cure.control <- function(tail_completion = c("zero", "exp", "tau-zero"),
                             tail_tau = Inf,
                             pmin = 1e-5,
                             ...)
{
    tail_comps <- c("zero", "exp", "tau-zero")
    tail_completion <- match.arg(tail_completion, choices = tail_comps)
    int_tail_completion <- match(tail_completion, tail_comps, nomatch = 1L)
    control0 <- do.call(intsurv.control, list(...))
    structure(c(control0,
                list(tail_completion = int_tail_completion,
                     tail_tau = tail_tau,
                     pmin = pmin)),
              class = c("cox_cure.control", "intsurv.control"))
}


##' @rdname intsurv.control
##'
##' @param nlambda A positive integer representing the number of lambda
##'     parameters.
##' @param lambda_min_ratio A positive number specifying the ratio between the
##'     smallest lambda in the solution path to the large enough lambda that
##'     would result in all zero estimates with the lasso penalty.
##' @param alpha A positive number between 0 and 1 representing the mixing
##'     parameter in the elastic net penalty.
##' @param lambda A numeric vector that consists of nonnegative values
##'     representing the sequence of the lambda parameters.
##' @param penalty_factor A numeric vector that consists of nonnegative penalty
##'     factors (or adaptive weights) for the L1-norm of the coefficient
##'     estimates.
##' @param varying_active A logical value.  If \code{TRUE} (by default), the
##'     underlying coordinate-descent algorithm will be iterated over varying
##'     active sets, which can usually improve the computational efficiency when
##'     the number of predictors is large.  Otherwise, an ordinary
##'     coordinate-descent will be performed.
##'
##' @export
cox_cure_net.control <- function(nlambda = 10,
                                 lambda_min_ratio = 1e-3,
                                 alpha = 1,
                                 lambda = NULL,
                                 penalty_factor = NULL,
                                 varying_active = TRUE,
                                 ...)
{
    control0 <- do.call(intsurv.control, list(...))
    structure(c(control0,
                list(nlambda = nlambda,
                     lambda_min_ratio = lambda_min_ratio,
                     alpha = alpha,
                     lambda = null2num0(lambda),
                     penalty_factor = null2num0(penalty_factor))),
              class = c("cox_cure_net.control", "intsurv.control"))
}
