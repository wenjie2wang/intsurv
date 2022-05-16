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


##' Regularized Cox Cure Rate Model
##'
##' For right-censored data, fit a regularized Cox cure rate model through
##' elastic-net penalty following Masud et al. (2018), and Zou and Hastie
##' (2005).  For right-censored data with uncertain event status,
##' fit the regularized Cox cure model proposed by Wang et al. (2020).  Without
##' regularization, the model reduces to the regular Cox cure rate model (Kuk
##' and Chen, 1992; Sy and Taylor, 2000)
##'
##' The model estimation procedure follows expectation maximization (EM)
##' algorithm.  Variable selection procedure through regularization by elastic
##' net penalty is developed based on cyclic coordinate descent and
##' majorization-minimization (MM) algorithm.
##'
##' @inheritParams cox_cure
##'
##' @param cv_nfolds A nonnegative integer representing the number of folds in
##'     cross-validation.
##'
##' @return
##'
##' \code{cox_cure_net} object for regular Cox cure rate model or
##' \code{cox_cure_net_mar} object for Cox cure rate model with uncertain
##' events.
##'
##' @references
##'
##' Kuk, A. Y. C., & Chen, C. (1992). A mixture model combining logistic
##' regression with proportional hazards regression. \emph{Biometrika}, 79(3),
##' 531--541.
##'
##' Masud, A., Tu, W., & Yu, Z. (2018). Variable selection for mixture and
##' promotion time cure rate models. \emph{Statistical methods in medical
##' research}, 27(7), 2185--2199.
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
##' Zou, H., & Hastie, T. (2005). Regularization and variable selection via the
##' elastic net. \emph{Journal of the Royal Statistical Society}: Series B
##' (Statistical Methodology), 67(2), 301--320.
##'
##' @seealso
##'
##' \code{\link{cox_cure}} for regular Cox cure rate model.
##'
##' @example inst/examples/cox_cure_net.R
##' @export
cox_cure_net <- function(surv_formula,
                         cure_formula,
                         time,
                         event,
                         data,
                         subset,
                         contrasts = NULL,
                         cv_nfolds = 0,
                         control = cox_cure.control(),
                         surv_control = cox_cure_net.control(),
                         cure_control = cox_cure_net.control(),
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
    call_list$cv_nfolds <- as.integer(cv_nfolds)
    ## cox model does not have an intercept
    surv_is_intercept <- colnames(call_list$surv_x) == "(Intercept)"
    surv_has_intercept <- any(surv_is_intercept)
    if ((ncol(surv_x) - as.integer(surv_has_intercept)) == 0L) {
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
    cure_standardize <- ! cure_only_intercept
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
        ## prepare additional mar/mcar layer
        call_list$mar_x <- model_list$mar$x
        mar_control <- cox_cure_net.control(
            lambda = 0,
            alpha = 0,
            start = numeric(0),
            offset = numeric(0),
            standardize = FALSE
        )
        names(mar_control) <- paste0("mar_", names(mar_control))
        mar_control$mar_intercept <- FALSE
        call_list <- c(call_list, mar_control)
        is_arg_valid <- names(call_list) %in%
            names(formals(rcpp_coxph_cure_mar_vs))
        call_list <- call_list[is_arg_valid]
        out <- do.call(rcpp_coxph_cure_mar_vs, call_list)
        ## add class
        class(out) <- "cox_cure_net_mar"
    } else {
        is_arg_valid <- names(call_list) %in% names(formals(rcpp_coxph_cure_vs))
        call_list <- call_list[is_arg_valid]
        out <- do.call(rcpp_coxph_cure_vs, call_list)
        ## add class
        class(out) <- "cox_cure_net"
    }
    ## add covariate names back
    if (! is.null(surv_var_names <- colnames(call_list$surv_x))) {
        colnames(out$surv_coef) <- surv_var_names
    } else {
        colnames(out$surv_coef) <- paste0("x", seq_along(out$surv_coef))
    }
    if (! is.null(cure_var_names <- colnames(call_list$cure_x))) {
        colnames(out$cure_coef) <- c(
        { if (call_list$cure_intercept) "(Intercept)" else NULL },
        cure_var_names
        )
    } else {
        colnames(out$cure_coef) <-
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


##' @rdname cox_cure_net
##'
##' @export
cox_cure_net.fit <- function(surv_x,
                             cure_x,
                             time,
                             event,
                             cure_intercept = TRUE,
                             cv_nfolds = 0,
                             control = cox_cure.control(),
                             surv_control = cox_cure_net.control(),
                             cure_control = cox_cure_net.control(),
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
    ## check time
    if (anyNA(time)) {
        stop("Found NA's in 'time'.")
    }
    ## check event
    if (all(event[! is.na(event)] < 1)) {
        stop("No event can be found.")
    }
    ## prepare a list to call the underlying model estimation function
    names(surv_control) <- paste0("surv_", names(surv_control))
    names(cure_control) <- paste0("cure_", names(cure_control))
    call_list <- c(control, surv_control, cure_control)
    call_list$time <- time
    call_list$event <- event
    call_list$surv_x <- surv_x
    call_list$cure_x <- cure_x
    call_list$bootstrap <- bootstrap
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
                 "the largest event time.")
        }
    }
    ## more checks if tail completion after a specified tau
    if (tail_completion == 3L) {
        if (tail_tau < max(time[! is.na(event) & event > 0])) {
            stop("The specified 'tail_tau' cannot be less than",
                 "the largest event time.")
        } else if (tail_tau > max(time)) {
            warning("The specified 'tail_tau' is greater than",
                    "the largest survival time.")
        }
    }
    ## call the routine
    if (anyNA(event)) {
        ## prepare additional mar/mcar layer
        call_list$mar_x <- matrix(1, nrow = length(time), ncol = 1)
        mar_control <- cox_cure_net.control(
            lambda = 0,
            offset = numeric(0),
            standardize = FALSE
        )
        names(mar_control) <- paste0("mar_", names(mar_control))
        mar_control$mar_intercept <- FALSE
        call_list <- c(call_list, mar_control)
        is_arg_valid <- names(call_list) %in%
            names(formals(rcpp_coxph_cure_mar_vs))
        call_list <- call_list[is_arg_valid]
        out <- do.call(rcpp_coxph_cure_mar_vs, call_list)
        ## add class
        class(out) <- "cox_cure_net_mar"
    } else {
        is_arg_valid <- names(call_list) %in% names(formals(rcpp_coxph_cure_vs))
        call_list <- call_list[is_arg_valid]
        out <- do.call(rcpp_coxph_cure_vs, call_list)
        ## add class
        class(out) <- "cox_cure_net"
    }
    ## add covariate names back
    if (! is.null(surv_var_names <- colnames(surv_x))) {
        colnames(out$surv_coef) <- surv_var_names
    } else {
        colnames(out$surv_coef) <- paste0("x", seq_len(ncol(out$surv_coef)))
    }
    if (! is.null(cure_var_names <- colnames(cure_x))) {
        colnames(out$cure_coef) <- c(
        {if (cure_intercept) "(Intercept)" else NULL},
        colnames(cure_x)
        )
    } else {
        colnames(out$cure_coef) <-
            if (cure_intercept) {
                c("(Intercept)",
                  paste0("z", seq_len(ncol(out$cure_coef) - 1L)))
            } else {
                paste0("z", seq_len(ncol(out$cure_coef)))
            }
    }
    if (control$save_call) {
        ## add function call
        out$call <- call0
    }
    ## return
    out
}
