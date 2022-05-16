##
## intsurv: Integrative Survival Models
## Copyright (C) 2017-2021  Wenjie Wang <wang@wwenjie.org>
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
##' 1992; Sy and Taylor, 2000) via an EM algorithm.
##'
##' @param surv_formula A formula object starting with \code{~} for the model
##'     formula in survival model part.  For Cox model, no intercept term is
##'     included even if an intercept is specified or implied in the model
##'     formula.  A model formula with an intercept term only is not allowed.
##' @param cure_formula A formula object starting with \code{~} for the model
##'     formula in incidence model part.  For logistic model, an intercept term
##'     is included by default and can be excluded by adding \code{+ 0} or
##'     \code{- 1} to the model formula.
##' @param time A numeric vector for the observed survival times.
##' @param event A numeric vector for the event indicators.
##' @param data An optional data frame, list, or environment that contains the
##'     model covariates and response variables (\code{time} and \code{event})
##'     If they are not found in data, the variables are taken from the
##'     environment of the specified formula, usually the environment from which
##'     this function is called.
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
##'     procedure will not run if \code{bootstrap = 0} by default.  If
##'     \code{bootstrap > 0}, the specified number of bootstrap samples will be
##'     used for estimating the standard errors.
##' @param control A \code{cox_cure.control} object that contains the control
##'     parameters.
##' @param surv_control,cure_control \code{intsurv.control} objects that contain
##'     the control parameters for the survival model part and logistic model
##'     part, respectively.
##' @param ... Other arguments passed to the control functions for backward
##'     compatibility.
##'
##' @return
##' A \code{cox_cure} object that contains the fitted regular Cox cure rate
##' model.
##'
##' @references
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
##' @seealso
##' \code{\link{cox_cure_net}} for fitting the regularized Cox cure rate model
##' with elastic-net penalty.
##'
##' @example inst/examples/cox_cure.R
##'
##' @export
cox_cure <- function(surv_formula,
                     cure_formula,
                     time,
                     event,
                     data,
                     subset,
                     contrasts = NULL,
                     bootstrap = 0,
                     control = cox_cure.control(),
                     surv_control = intsurv.control(),
                     cure_control = intsurv.control(),
                     ...)
{
    ## controls
    if (! inherits(control, "cox_cure.control")) {
        control <- do.call(cox_cure.control, control)
    }
    if (! inherits(surv_control, "intsurv.control")) {
        surv_control <- do.call(intsurv.control, surv_control)
    }
    if (! inherits(cure_control, "intsurv.control")) {
        cure_control <- do.call(intsurv.control, cure_control)
    }
    if (control$save_call) {
        ## record function call
        call0 <- match.call()
    }
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
    ## prepare a list to call the underlying model estimation function
    names(surv_control) <- paste0("surv_", names(surv_control))
    names(cure_control) <- paste0("cure_", names(cure_control))
    call_list <- c(control, surv_control, cure_control)
    ## get design matrix and responses
    call_list$time <- model_list$surv$time
    call_list$event <- model_list$surv$event
    call_list$surv_x <- model_list$surv$x
    call_list$cure_x <- model_list$cure$x
    call_listsurv_offset <- model_list$surv$offset
    call_list$cure_offset <- model_list$cure$offset
    call_list$bootstrap <- as.integer(bootstrap)
    ## cox model does not have an intercept
    surv_is_intercept <- colnames(call_list$surv_x) == "(Intercept)"
    surv_has_intercept <- any(surv_is_intercept)
    if ((ncol(call_list$surv_x) - as.integer(surv_has_intercept)) == 0L) {
        stop("No covariate is specified in 'formula'.")
    }
    ## remove the possible intercept term
    if (surv_has_intercept) {
        call_list$surv_x <- call_list$surv_x[, which(! surv_is_intercept),
                                             drop = FALSE]
    }
    ## logistic model can have an intercept
    cure_is_intercept <- colnames(call_list$cure_x) == "(Intercept)"
    cure_has_intercept <- any(cure_is_intercept)
    cure_only_intercept <- all(cure_is_intercept)
    call_list$cure_standardize <- ! cure_only_intercept
    ## remove the possible intercept term
    if (cure_has_intercept) {
        call_list$cure_x <- call_list$cure_x[, which(! cure_is_intercept),
                                             drop = FALSE]
        call_list$cure_intercept <- TRUE
    } else {
        call_list$cure_intercept <- FALSE
    }
    ## check event
    if (all(call_list$event[! is.na(call_list$event)] < 1)) {
        stop("No event can be found.")
    }
    ## start values
    call_list$surv_start <- null2num0(call_list$surv_start)
    call_list$cure_start <- null2num0(call_list$cure_start)
    ## offset terms
    call_list$surv_offset <- null2num0(call_list$surv_offset)
    call_list$cure_offset <- null2num0(call_list$cure_offset)
    ## more checks if tail completion after a specified tau
    if (call_list$tail_completion == 3L) {
        is_tau_small <- with(call_list,
                             tail_tau < max(time[! is.na(event) & event > 0]))
        if (is_tau_small) {
            stop("The specified 'tail_tau' cannot be less than",
                 " the largest event time.")
        }
    }
    ## call the routine
    if (anyNA(call_list$event)) {
        ## prepare additional mar
        call_list$mar_x <- model_list$mar$x
        mar_control <- intsurv.control(
            start = numeric(0),
            offset = numeric(0),
            standardize = FALSE
        )
        names(mar_control) <- paste0("mar_", names(mar_control))
        mar_control$mar_intercept <- FALSE
        call_list <- c(call_list, mar_control)
        is_arg_valid <- names(call_list) %in%
            names(formals(rcpp_coxph_cure_mar))
        call_list <- call_list[is_arg_valid]
        out <- do.call(rcpp_coxph_cure_mar, call_list)
        ## add class
        class(out) <- "cox_cure_mar"
    } else {
        is_arg_valid <- names(call_list) %in% names(formals(rcpp_coxph_cure))
        call_list <- call_list[is_arg_valid]
        out <- do.call(rcpp_coxph_cure, call_list)
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
    if (! is.null(surv_var_names <- colnames(call_list$surv_x))) {
        names(out$surv_coef) <- surv_var_names
    } else {
        names(out$surv_coef) <- paste0("x", seq_along(out$surv_coef))
    }
    if (! is.null(cure_var_names <- colnames(call_list$cure_x))) {
        names(out$cure_coef) <- c(
        { if (call_list$cure_intercept) "(Intercept)" else NULL },
        cure_var_names
        )
    } else {
        names(out$cure_coef) <-
            if (call_list$cure_intercept) {
                c("(Intercept)",
                  paste0("z", seq_along(out$cure_coef[- 1L])))
            } else {
                paste0("z", seq_along(out$cure_coef))
            }
    }
    if (control$save_call) {
        ## add function call
        out$call <- call0
    }
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
cox_cure.fit <- function(surv_x,
                         cure_x,
                         time,
                         event,
                         cure_intercept = TRUE,
                         bootstrap = 0,
                         control = cox_cure.control(),
                         surv_control = intsurv.control(),
                         cure_control = intsurv.control(),
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
                    "the largest survival time")
        }
    }

    ## call the routine
    if (anyNA(event)) {
        out <- coxph_cure_mar(
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
            em_rel_tol = em_rel_tol,
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
        class(out) <- "cox_cure_mcar"
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
            em_rel_tol = em_rel_tol,
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
