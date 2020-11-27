##
## intsurv: Integrative Survival Models
## Copyright (C) 2017-2020  Wenjie Wang <wang@wwenjie.org>
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


##' Cox Cure Rate Model
##'
##' For right-censored data, fit a regular Cox cure rate model (Kuk and Chen,
##' 1992; Sy and Taylor, 2000) via an EM algorithm.  For right-censored data
##' with uncertain event status, fit the Cox cure model proposed by Wang et
##' al. (2020).
##'
##' @param surv_formula A formula object starting with \code{~} for the model
##'     formula in survival model part.  For Cox model, no intercept term is
##'     included even if an intercept is specified or implied in the model
##'     formula.  A model formula with an intercept term only is not allowed.
##' @param cure_formula A formula object starting with \code{~} for the model
##'     formula in cure rate model part.  For logistic model, an intercept term
##'     is included by default and can be excluded by adding \code{+ 0} or
##'     \code{- 1} to the model formula.
##' @param time A numeric vector for the observed survival times.
##' @param event A numeric vector for the event indicators.  \code{NA}'s are
##'     allowed and represent uncertain event status.
##' @param data An optional data frame, list, or environment that contains the
##'     covariates and response variables (\code{time} and \code{event}),
##'     included in the model. If they are not found in data, the variables are
##'     taken from the environment of the specified formula, usually the
##'     environment from which this function is called.
##' @param subset An optional logical vector specifying a subset of observations
##'     to be used in the fitting process.
##' @param contrasts An optional list, whose entries are values (numeric
##'     matrices or character strings naming functions) to be used as
##'     replacement values for the contrasts replacement function and whose
##'     names are the names of columns of data containing factors.  See
##'     \code{contrasts.arg} of \code{\link[stats]{model.matrix.default}} for
##'     details.
##' @param bootstrap An integer representing the number of bootstrap samples for
##'     estimating standard errors of the coefficient estimates.  The bootstrap
##'     procedure will not run with \code{bootstrap = 0} by default.  If
##'     \code{bootstrap > 0}, the specified number of bootstrap samples will be
##'     used in estimating standard errors.
##' @param firth A logical value indicating whether to use Firth's
##'     bias-reduction method (Firth, 1993) in the logistic model component.
##'     The default value is \code{FALSE} for fitting a regular logistic model.
##'     Notice that this argument is experimental and only available for regular
##'     Cox cure rate model currently.
##' @param surv_start An optional numeric vector representing the starting
##'     values for the Cox model component.  If \code{NULL} is specified, the
##'     starting values will be obtained from fitting a regular Cox model to
##'     events only.
##' @param cure_start An optional numeric vector representing the starting
##'     values for the logistic model component.  If \code{NULL} is specified,
##'     the starting values will be obtained from fitting a regular logistic
##'     model to the non-missing event indicators.
##' @param surv_offset An optional numeric vector representing the offset term
##'     in the Cox model compoent.
##' @param cure_offset An optional numeric vector representing the offset term
##'     in the logistic model compoent.
##' @param em_max_iter A positive integer specifying the maximum iteration
##'     number of the EM algorithm.  The default value is \code{200}.
##' @param em_rel_tol A positive number specifying the tolerance that determines
##'     the convergence of the EM algorithm in terms of the convergence of the
##'     covariate coefficient estimates.  The tolerance is compared with the
##'     relative change between estimates from two consecutive iterations, which
##'     is measured by ratio of the L1-norm of their difference to the sum of
##'     their L1-norm.  The default value is \code{1e-5}.
##' @param surv_max_iter A positive integer specifying the maximum iteration
##'     number of the M-step routine related to the survival model component.
##'     The default value is \code{200}.
##' @param surv_rel_tol A positive number specifying the tolerance that
##'     determines the convergence of the M-step related to the survival model
##'     component in terms of the convergence of the covariate coefficient
##'     estimates.  The tolerance is compared with the relative change between
##'     estimates from two consecutive iterations, which is measured by ratio of
##'     the L1-norm of their difference to the sum of their L1-norm.  The
##'     default value is \code{1e-5}.
##' @param cure_max_iter A positive integer specifying the maximum iteration
##'     number of the M-step routine related to the cure rate model component.
##'     The default value is \code{200}.
##' @param cure_rel_tol A positive number specifying the tolerance that
##'     determines the convergence of the M-step related to the cure rate model
##'     component in terms of the convergence of the covariate coefficient
##'     estimates.  The tolerance is compared with the relative change between
##'     estimates from two consecutive iterations, which is measured by ratio of
##'     the L1-norm of their difference to the sum of their L1-norm.  The
##'     default value is \code{1e-5}.
##' @param tail_completion A character string specifying the tail completion
##'     method for conditional survival function.  The available methods are
##'     \code{"zero"} for zero-tail completion after the largest event times (Sy
##'     and Taylor, 2000), \code{"exp"} for exponential-tail completion (Peng,
##'     2003), and \code{"zero-tau"} for zero-tail completion after a specified
##'     \code{tail_tau}.  The default method is the zero-tail completion
##'     proposed by Sy and Taylor (2000).
##' @param tail_tau A numeric number specifying the time of zero-tail
##'     completion.  It will be used only if \code{tail_completion =
##'     "zero-tau"}.  A reasonable choice must be a time point between the
##'     largest event time and the largest survival time.
##' @param pmin A numeric number specifying the minimum value of probabilities
##'     for sake of numerical stability.  The default value is \code{1e-5}.
##' @param early_stop A logical value specifying whether to stop the iteration
##'     once the negative log-likelihood unexpectedly increases, which may
##'     suggest convergence on likelihood, or indicate numerical issues or
##'     implementation bugs.  The default value is \code{TRUE}.
##' @param verbose A logical value.  If \code{TRUE}, a verbose information will
##'     be given along iterations for tracing the convergence.  The default
##'     value is \code{FALSE}.
##' @param ... Other arguments for future usage.  A warning will be thrown if
##'     any invalid argument is specified.
##'
##' @return
##'
##' \code{cox_cure} object for regular Cox cure rate model or
##' \code{cox_cure_uncer} object for Cox cure rate model with uncertain events.
##'
##' @references
##'
##' Firth, D. (1993). Bias reduction of maximum likelihood
##' estimates. \emph{Biometrika}, 80(1), 27--38.
##'
##' Kuk, A. Y. C., & Chen, C. (1992). A mixture model combining logistic
##' regression with proportional hazards regression. \emph{Biometrika}, 79(3),
##' 531--541.
##'
##' Peng, Y. (2003). Estimating baseline distribution in proportional hazards
##' cure models. \emph{Computational Statistics & Data Analysis}, 42(1-2),
##' 187--201.
##'
##' Sy, J. P., & Taylor, J. M. (2000). Estimation in a Cox proportional hazards
##' cure model. \emph{Biometrics}, 56(1), 227--236.
##'
##' Wang, W., Luo, C., Aseltine, R. H., Wang, F., Yan, J., & Chen, K. (2020).
##' Suicide Risk Modeling with Uncertain Diagnostic Records. \emph{arXiv
##' preprint arXiv:2009.02597}.
##'
##' @seealso
##'
##' \code{\link{cox_cure_net}} for regularized Cox cure rate model with
##' elastic-net penalty.
##'
##' @example inst/examples/cox_cure.R
##'
##' @export
cox_cure <- function(surv_formula,
                     cure_formula,
                     time, event,
                     data, subset,
                     contrasts = NULL,
                     bootstrap = 0,
                     firth = FALSE,
                     surv_start = NULL,
                     cure_start = NULL,
                     surv_offset = NULL,
                     cure_offset = NULL,
                     em_max_iter = 200,
                     em_rel_tol = 1e-5,
                     surv_max_iter = 30,
                     surv_rel_tol = 1e-5,
                     cure_max_iter = 30,
                     cure_rel_tol = 1e-5,
                     tail_completion = c("zero", "exp", "zero-tau"),
                     tail_tau = NULL,
                     pmin = 1e-5,
                     early_stop = TRUE,
                     verbose = FALSE,
                     ...)
{
    ## warning on `...`
    warn_dots(...)

    ## record function call
    call0 <- match.call()

    ## prepare to call function prep_cure_model
    this_call <- match.call(expand.dots = FALSE)
    ## time is also a function name.  rename to avoid potential issues.
    names(this_call)[which(names(this_call) == "time")] <- "obs_time"
    names(this_call)[which(names(this_call) == "event")] <- "obs_event"
    this_call$eval_env <- parent.frame()
    matched_call <- match(names(formals(prep_cure_model)),
                          names(this_call), nomatch = 0L)
    model_call <- this_call[c(1L, matched_call)]
    model_call[[1L]] <- quote(prep_cure_model)
    model_list <- eval(model_call)

    ## get design matrix and responses
    obs_time <- model_list$surv$obs_time
    obs_event <- model_list$surv$obs_event
    surv_x <- model_list$surv$xMat
    cure_x <- model_list$cure$xMat

    ## cox model does not have an intercept
    surv_is_intercept <- colnames(surv_x) == "(Intercept)"
    surv_has_intercept <- any(surv_is_intercept)
    if ((ncol(surv_x) - as.integer(surv_has_intercept)) == 0L) {
        stop("No covariate is specified in 'formula'.")
    }
    ## remove the possible intercept term
    if (surv_has_intercept) {
        surv_x <- surv_x[, which(! surv_is_intercept), drop = FALSE]
    }

    ## logistic model can have an intercept
    cure_is_intercept <- colnames(cure_x) == "(Intercept)"
    cure_has_intercept <- any(cure_is_intercept)
    cure_only_intercept <- all(cure_is_intercept)
    cure_standardize <- ! cure_only_intercept
    ## remove the possible intercept term
    if (cure_has_intercept) {
        cure_x <- cure_x[, which(! cure_is_intercept), drop = FALSE]
        cure_intercept <- TRUE
    } else {
        cure_intercept <- FALSE
    }

    ## check time
    if (anyNA(obs_time)) {
        stop("Found NA's in 'time'.")
    }
    ## check event
    if (all(obs_event[! is.na(obs_event)] < 1)) {
        stop("No event can be found.")
    }
    ## starting values
    if (is.null(surv_start)) {
        surv_start <- 0
    } else if (length(surv_start) != ncol(surv_x)) {
        stop("The length of 'surv_start' is inappropriate.")
    }
    if (is.null(cure_start)) {
        cure_start <- 0
    } else if (length(cure_start) != ncol(cure_x) +
               as.integer(cure_intercept)) {
        stop("The length of 'cure_start' is inappropriate.")
    }
    ## offset terms
    if (is.null(surv_offset)) {
        surv_offset <- rep(0, nrow(surv_x))
    } else if (length(surv_offset) != nrow(surv_x)) {
        stop("The length of 'surv_offset' is inappropriate.")
    }
    if (is.null(cure_offset)) {
        cure_offset <- rep(0, nrow(cure_x))
    } else if (length(cure_offset) != nrow(cure_x)) {
        stop("The length of 'cure_start' is inappropriate.")
    }
    ## on tail completion
    all_tails <- c("zero", "exp", "zero-tau")
    tail_completion <- match(match.arg(tail_completion, all_tails),
                             all_tails)
    ## prepare tail_tau
    if (is.null(tail_tau)) {
        tail_tau <- - 1
    }
    ## more checks if tail completion after a specified tau
    if (tail_completion == 3L) {
        if (tail_tau < max(obs_time[! is.na(obs_event) & obs_event > 0])) {
            stop("The specified 'tail_tau' cannot be less than",
                 "the largest event time.")
        } else if (tail_tau > max(obs_time)) {
            warning("The specified 'tail_tau' is greater than",
                    "the largest survival time, which ")
        }
    }

    ## call the routine
    if (anyNA(obs_event)) {
        out <- coxph_cure_uncer(
            time = obs_time,
            event = obs_event,
            cox_x = surv_x,
            cure_x = cure_x,
            cure_intercept = cure_intercept,
            bootstrap = bootstrap,
            cox_start = surv_start,
            cure_start = cure_start,
            cox_offset = surv_offset,
            cure_offset = cure_offset,
            cox_standardize = TRUE,
            cure_standardize = cure_standardize,
            em_max_iter = em_max_iter,
            cox_mstep_max_iter = surv_max_iter,
            cox_mstep_rel_tol = surv_rel_tol,
            cure_mstep_max_iter = cure_max_iter,
            cure_mstep_rel_tol = cure_rel_tol,
            tail_completion = tail_completion,
            tail_tau = tail_tau,
            pmin = pmin,
            early_stop = early_stop,
            verbose = verbose
        )
        ## add class
        class(out) <- "cox_cure_uncer"
    } else {
        out <- rcpp_coxph_cure(
            time = obs_time,
            event = obs_event,
            cox_x = surv_x,
            cure_x = cure_x,
            cure_intercept = cure_intercept,
            bootstrap = bootstrap,
            firth = firth,
            cox_start = surv_start,
            cure_start = cure_start,
            cox_offset = surv_offset,
            cure_offset = cure_offset,
            cox_standardize = TRUE,
            cure_standardize = cure_standardize,
            em_max_iter = em_max_iter,
            cox_mstep_max_iter = surv_max_iter,
            cox_mstep_rel_tol = surv_rel_tol,
            cure_mstep_max_iter = cure_max_iter,
            cure_mstep_rel_tol = cure_rel_tol,
            tail_completion = tail_completion,
            tail_tau = tail_tau,
            pmin = pmin,
            early_stop = early_stop,
            verbose = verbose
        )
        ## add class
        class(out) <- "cox_cure"
    }
    ## add bootstrap se if available
    if (bootstrap <= 0) {
        out$bootstrap$surv_coef_mat <- out$bootstrap$cure_coef_mat <- NULL
    } else {
        out$bootstrap$surv_coef_se <- apply(out$bootstrap$surv_coef_mat,
                                            2L, se_interQ)
        out$bootstrap$cure_coef_se <- apply(out$bootstrap$cure_coef_mat,
                                            2L, se_interQ)
    }
    ## add covariate names back
    if (! is.null(surv_var_names <- colnames(surv_x))) {
        names(out$surv_coef) <- surv_var_names
    } else {
        names(out$surv_coef) <- paste0("x", seq_along(out$surv_coef))
    }
    if (! is.null(cure_var_names <- colnames(cure_x))) {
        names(out$cure_coef) <- c({if (cure_intercept) "(Intercept)" else NULL},
                                  colnames(cure_x))
    } else {
        names(out$cure_coef) <-
            if (cure_intercept) {
                c("(Intercept)",
                  paste0("z", seq_along(out$cure_coef[- 1L])))
            } else {
                paste0("z", seq_along(out$cure_coef))
            }
    }
    ## add function call
    out$call <- call0
    ## return
    out
}


##' @rdname cox_cure
##'
##' @param surv_x A numeric matrix for the design matrix of the survival model
##'     component.
##' @param cure_x A numeric matrix for the design matrix of the cure rate model
##'     component.  The design matrix should exclude an intercept term unless we
##'     want to fit a model only including the intercept term.  In that case, we
##'     need further set \code{cure_intercept = FALSE} to not standardize the
##'     intercept term.
##' @param cure_intercept A logical value specifying whether to add an intercept
##'     term to the cure rate model component.  If \code{TRUE} by default, an
##'     intercept term is included.
##' @param surv_standardize A logical value specifying whether to standardize
##'     the covariates for the survival model part.  If \code{FALSE}, the
##'     covariates will be standardized internally to have mean zero and
##'     standard deviation one.
##' @param cure_standardize A logical value specifying whether to standardize
##'     the covariates for the cure rate model part.  If \code{TRUE} (by
##'     default), the covariates will be standardized internally to have mean
##'     zero and standard deviation one.
##'
##' @export
cox_cure.fit <- function(surv_x, cure_x,
                         time, event,
                         cure_intercept = TRUE,
                         bootstrap = 0,
                         firth = FALSE,
                         surv_start = NULL,
                         cure_start = NULL,
                         surv_offset = NULL,
                         cure_offset = NULL,
                         surv_standardize = TRUE,
                         cure_standardize = TRUE,
                         em_max_iter = 200,
                         em_rel_tol = 1e-5,
                         surv_max_iter = 30,
                         surv_rel_tol = 1e-5,
                         cure_max_iter = 30,
                         cure_rel_tol = 1e-5,
                         tail_completion = c("zero", "exp", "zero-tau"),
                         tail_tau = NULL,
                         pmin = 1e-5,
                         early_stop = TRUE,
                         verbose = FALSE,
                         ...)
{
    ## warning on `...`
    warn_dots(...)

    ## record function call
    this_call <- match.call()

    ## check time
    if (anyNA(time)) {
        stop("Found NA's in 'time'.")
    }
    ## check event
    if (all(event[! is.na(event)] < 1)) {
        stop("No event can be found.")
    }
    ## starting values
    if (is.null(surv_start)) {
        surv_start <- 0
    } else if (length(surv_start) != ncol(surv_x)) {
        stop("The length of 'surv_start' is inappropriate.")
    }
    if (is.null(cure_start)) {
        cure_start <- 0
    } else if (length(cure_start) != ncol(cure_x) +
               as.integer(cure_intercept)) {
        stop("The length of 'cure_start' is inappropriate.")
    }
    ## offset terms
    if (is.null(surv_offset)) {
        surv_offset <- rep(0, nrow(surv_x))
    } else if (length(surv_offset) != nrow(surv_x)) {
        stop("The length of 'surv_offset' is inappropriate.")
    }
    if (is.null(cure_offset)) {
        cure_offset <- rep(0, nrow(cure_x))
    } else if (length(cure_offset) != nrow(cure_x)) {
        stop("The length of 'cure_start' is inappropriate.")
    }
    ## on tail completion
    all_tails <- c("zero", "exp", "zero-tau")
    tail_completion <- match(match.arg(tail_completion, all_tails),
                             all_tails)
    ## prepare tail_tau
    if (is.null(tail_tau)) {
        tail_tau <- - 1
    }
    ## more checks if tail completion after a specified tau
    if (tail_completion == 3L) {
        if (tail_tau < max(time[! is.na(event) & event > 0])) {
            stop("The specified 'tail_tau' cannot be less than",
                 "the largest event time.")
        } else if (tail_tau > max(time)) {
            warning("The specified 'tail_tau' is greater than",
                    "the largest survival time, which ")
        }
    }

    ## call the routine
    if (anyNA(event)) {
        out <- coxph_cure_uncer(
            time = time,
            event = event,
            cox_x = surv_x,
            cure_x = cure_x,
            cure_intercept = cure_intercept,
            bootstrap = bootstrap,
            cox_start = surv_start,
            cure_start = cure_start,
            cox_offset = surv_offset,
            cure_offset = cure_offset,
            cox_standardize = surv_standardize,
            cure_standardize = cure_standardize,
            em_max_iter = em_max_iter,
            cox_mstep_max_iter = surv_max_iter,
            cox_mstep_rel_tol = surv_rel_tol,
            cure_mstep_max_iter = cure_max_iter,
            cure_mstep_rel_tol = cure_rel_tol,
            tail_completion = tail_completion,
            tail_tau = tail_tau,
            pmin = pmin,
            early_stop = early_stop,
            verbose = verbose
        )
        ## add class
        class(out) <- "cox_cure_uncer"
    } else {
        out <- rcpp_coxph_cure(
            time = time,
            event = event,
            cox_x = surv_x,
            cure_x = cure_x,
            cure_intercept = cure_intercept,
            bootstrap = bootstrap,
            firth = firth,
            cox_start = surv_start,
            cure_start = cure_start,
            cox_offset = surv_offset,
            cure_offset = cure_offset,
            cox_standardize = surv_standardize,
            cure_standardize = cure_standardize,
            em_max_iter = em_max_iter,
            cox_mstep_max_iter = surv_max_iter,
            cox_mstep_rel_tol = surv_rel_tol,
            cure_mstep_max_iter = cure_max_iter,
            cure_mstep_rel_tol = cure_rel_tol,
            tail_completion = tail_completion,
            tail_tau = tail_tau,
            pmin = pmin,
            early_stop = early_stop,
            verbose = verbose
        )
        ## add class
        class(out) <- "cox_cure"
    }
    ## add bootstrap se if available
    if (bootstrap <= 0) {
        out$bootstrap$surv_coef_mat <- out$bootstrap$cure_coef_mat <- NULL
    } else {
        out$bootstrap$surv_coef_se <- apply(out$bootstrap$surv_coef_mat,
                                            2L, se_interQ)
        out$bootstrap$cure_coef_se <- apply(out$bootstrap$cure_coef_mat,
                                            2L, se_interQ)
    }
    ## add covariate names back
    if (! is.null(surv_var_names <- colnames(surv_x))) {
        names(out$surv_coef) <- surv_var_names
    } else {
        names(out$surv_coef) <- paste0("x", seq_along(out$surv_coef))
    }
    if (! is.null(cure_var_names <- colnames(cure_x))) {
        names(out$cure_coef) <- c({if (cure_intercept) "(Intercept)" else NULL},
                                  colnames(cure_x))
    } else {
        names(out$cure_coef) <-
            if (cure_intercept) {
                c("(Intercept)",
                  paste0("z", seq_along(out$cure_coef[- 1L])))
            } else {
                paste0("z", seq_along(out$cure_coef))
            }
    }
    ## add function call
    out$call <- this_call
    ## return
    out
}
