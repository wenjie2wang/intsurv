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


##' Cox Cure Rate Model
##'
##' For right-censored data, the function \code{cox_cure()} trains a Cox cure
##' rate model (Kuk and Chen, 1992; Sy and Taylor, 2000) via an expectation
##' maximization (EM) algorithm; For right-censored data with missing/uncertain
##' event/censoring indicators, the function fits the Cox cure rate model
##' proposed by Wang et al. (2023).
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
##' @param event A numeric vector for the event indicators, where \code{NA}'s
##'     are allowed and represent uncertain event indicators.
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
##' @param surv_mstep,cure_mstep A named list passed to \code{cox_cure.mstep()}
##'     specifying the control parameters for the corresponding M-steps.
##' @param control A \code{cox_cure.control} object that contains the control
##'     parameters.
##' @param ... Other arguments passed to the control functions for backward
##'     compatibility.
##'
##' @return A \code{cox_cure} object that contains the fitted ordinary Cox cure
##'     rate model if none of the event indicators is \code{NA}.  For
##'     right-censored data with uncertain/missing event indicators, a
##'     \code{cox_cure_uncer} object is returned.
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
##' Wang, W., Luo, C., Aseltine, R. H., Wang, F., Yan, J., & Chen,
##' K. (2023). Survival Modeling of Suicide Risk with Rare and Uncertain
##' Diagnoses. \emph{Statistics in Biosciences}, 17(1), 1--27.
##'
##' @example inst/examples/ex-cox_cure.R
##'
##' @export
cox_cure <- function(surv_formula,
                     cure_formula,
                     time,
                     event,
                     data,
                     subset,
                     contrasts = NULL,
                     bootstrap = 0L,
                     surv_mstep = cox_cure.mstep(),
                     cure_mstep = cox_cure.mstep(),
                     control = cox_cure.control(),
                     ...)
{
    ## for backward compatibility
    dot_list <- list(...)
    update_ctrl <- function(ctrl, old_args, prefix = "surv_") {
        for (old_a in old_args) {
            old_a2 <- paste0(prefix, old_a)
            if (!is.null(dot_list[[old_a2]])) {
                warning(wrapMessages(sprintf("The argument '%s'", old_a2),
                                     "is deprecated for cox_cure()."),
                        call. = FALSE)
                ctrl[[old_a]] <- dot_list[[old_a2]]
            }
        }
        ctrl
    }
    surv_mstep <- update_ctrl(
        surv_mstep,
        old_args = c("start", "offset", "max_iter", "rel_tol", "standardize"),
        prefix = "surv_"
    )
    cure_mstep <- update_ctrl(
        cure_mstep,
        old_args = c("start", "offset", "max_iter", "rel_tol", "standardize"),
        prefix = "cure_"
    )
    control <- update_ctrl(
        control,
        old_args = c("max_iter", "rel_tol"),
        prefix = "em_"
    )
    control <- update_ctrl(
        control,
        old_args = c("tail_completion", "tail_tau", "pmin", "verbose"),
        prefix = ""
    )
    ## for others
    invalid_args <- c("firth", "early_stop")
    for (a in invalid_args) {
        if (!is.null(dot_list[[a]])) {
            warning(wrapMessages(sprintf(
                "The argument '%s' was experimental", a
            ), "and is no longer in use."))
        }
    }

    ## refresh controls
    control <- do.call(cox_cure.control, control)
    surv_mstep <- do.call(cox_cure.mstep, surv_mstep)
    cure_mstep <- do.call(cox_cure.mstep, cure_mstep)
    out_control <- list(
        surv_mstep = surv_mstep,
        cure_mstep = cure_mstep,
        control = control
    )
    ## prepare a list to call the underlying model estimation function
    names(surv_mstep) <- paste0("surv_", names(surv_mstep))
    names(cure_mstep) <- paste0("cure_", names(cure_mstep))
    call_list <- c(control, surv_mstep, cure_mstep)

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
    ## get design matrix and responses
    call_list$time <- model_list$surv$time
    call_list$event <- model_list$surv$event
    call_list$surv_x <- model_list$surv$x
    call_list$cure_x <- model_list$cure$x
    if (! is.null(model_list$surv$offset)) {
        call_list$surv_offset <- model_list$surv$offset
    }
    if (! is.null(model_list$cure$offset)) {
        call_list$cure_offset <- model_list$cure$offset
    }
    call_list$bootstrap <- as.integer(bootstrap)
    ## cox model does not have an intercept
    surv_is_intercept <- colnames(call_list$surv_x) == "(Intercept)"
    surv_has_intercept <- any(surv_is_intercept)
    if ((ncol(call_list$surv_x) - as.integer(surv_has_intercept)) == 0L) {
        stop("No covariate is specified in 'surv_formula'.")
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
    has_uncer <- anyNA(call_list$event)
    if (has_uncer) {
        is_arg_valid <- names(call_list) %in%
            names(formals(rcpp_coxph_cure_mar))
        call_list <- call_list[is_arg_valid]
        out <- do.call(rcpp_coxph_cure_mar, call_list)
    } else {
        is_arg_valid <- names(call_list) %in% names(formals(rcpp_coxph_cure))
        call_list <- call_list[is_arg_valid]
        out <- do.call(rcpp_coxph_cure, call_list)
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
    ## add controls
    out <- c(out, out_control)
    ## return
    if (has_uncer) {
        structure(out, class = "cox_cure_uncer")
    } else {
        structure(out, class = "cox_cure")
    }
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
##'
##' @export
cox_cure.fit <- function(surv_x,
                         cure_x,
                         time,
                         event,
                         cure_intercept = TRUE,
                         bootstrap = 0L,
                         surv_mstep = cox_cure.mstep(),
                         cure_mstep = cox_cure.mstep(),
                         control = cox_cure.control(),
                         ...)
{
    ## for backward compatibility
    dot_list <- list(...)
    update_ctrl <- function(ctrl, old_args, prefix = "surv_") {
        for (old_a in old_args) {
            old_a2 <- paste0(prefix, old_a)
            if (!is.null(dot_list[[old_a2]])) {
                warning(wrapMessages(sprintf("The argument '%s'", old_a2),
                                     "is deprecated for cox_cure.fit()."),
                        call. = FALSE)
                ctrl[[old_a]] <- dot_list[[old_a2]]
            }
        }
        ctrl
    }
    surv_mstep <- update_ctrl(
        surv_mstep,
        old_args = c("start", "offset", "max_iter", "rel_tol", "standardize"),
        prefix = "surv_"
    )
    cure_mstep <- update_ctrl(
        cure_mstep,
        old_args = c("start", "offset", "max_iter", "rel_tol", "standardize"),
        prefix = "cure_"
    )
    control <- update_ctrl(
        control,
        old_args = c("max_iter", "rel_tol"),
        prefix = "em_"
    )
    control <- update_ctrl(
        control,
        old_args = c("tail_completion", "tail_tau", "pmin", "verbose"),
        prefix = ""
    )

    ## for others
    invalid_args <- c("firth", "early_stop")
    for (a in invalid_args) {
        if (!is.null(dot_list[[a]])) {
            warning(wrapMessages(sprintf(
                "The argument '%s' was experimental", a
            ), "and is no longer in use."))
        }
    }

    ## refresh controls
    control <- do.call(cox_cure.control, control)
    surv_mstep <- do.call(cox_cure.mstep, surv_mstep)
    cure_mstep <- do.call(cox_cure.mstep, cure_mstep)
    out_control <- list(
        surv_mstep = surv_mstep,
        cure_mstep = cure_mstep,
        control = control
    )
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
    names(surv_mstep) <- paste0("surv_", names(surv_mstep))
    names(cure_mstep) <- paste0("cure_", names(cure_mstep))
    call_list <- c(control, surv_mstep, cure_mstep)
    call_list$time <- time
    call_list$event <- event
    call_list$surv_x <- surv_x
    call_list$cure_x <- cure_x
    call_list$bootstrap <- bootstrap
    call_list$cure_intercept <- cure_intercept
    ## more checks if tail completion after a specified tau
    if (call_list$tail_completion == 3L) {
        is_tau_small <- with(call_list,
                             tail_tau < max(time[! is.na(event) & event > 0]))
        if (is_tau_small) {
            stop("The specified 'tail_tau' cannot be less than",
                 "the largest event time.")
        }
    }
    ## call the routine
    has_uncer <- anyNA(event)
    if (has_uncer) {
        is_arg_valid <- names(call_list) %in%
            names(formals(rcpp_coxph_cure_mar))
        call_list <- call_list[is_arg_valid]
        out <- do.call(rcpp_coxph_cure_mar, call_list)
    } else {
        is_arg_valid <- names(call_list) %in% names(formals(rcpp_coxph_cure))
        call_list <- call_list[is_arg_valid]
        out <- do.call(rcpp_coxph_cure, call_list)
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
                                  cure_var_names)
    } else {
        names(out$cure_coef) <-
            if (cure_intercept) {
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
    ## add controls
    out <- c(out, out_control)
    ## return
    if (has_uncer) {
        structure(out, class = "cox_cure_uncer")
    } else {
        structure(out, class = "cox_cure")
    }
}


##' @rdname cox_cure
##'
##' @param maxit A positive integer specifying the maximum iteration number.
##'     The default value is \code{1000}.
##' @param epsilon A positive number specifying the tolerance that determines
##'     the convergence of the coefficient estimates.  The tolerance is compared
##'     with the relative change between estimates from two consecutive
##'     iterations, which is measured by the ratio of the L1-norm of their
##'     difference to the sum of their L1-norms plus one.
##' @param verbose A nonnegative integer for verbose outputs, which is mainly
##'     useful for debugging.
##' @param tail_completion A character string specifying the tail completion
##'     method for conditional survival function.  The available methods are
##'     \code{"zero"} for zero-tail completion after the largest event times (Sy
##'     and Taylor, 2000), \code{"exp"} for exponential-tail completion (Peng,
##'     2003), and \code{"tau-zero"} for zero-tail completion after a specified
##'     \code{tail_tau}.  The default method is the zero-tail completion
##'     proposed by Sy and Taylor (2000).
##' @param tail_tau A numeric number specifying the time of zero-tail
##'     completion.  It will be used only if \code{tail_completion =
##'     "tau-zero"}.  A reasonable choice must be a time point between the
##'     largest event time and the largest survival time.
##' @param pmin A positive number specifying the minimum value of probabilities
##'     for numerical stability.  The default value is \code{1e-5}.
##' @param save_call A logical value indicating if the function call should be
##'     saved.  For large datasets, saving the function call would increase the
##'     size of the returned object dramatically.  We may want to set
##'     \code{save_call = FALSE} if the original function call is not needed.
##'
##' @export
cox_cure.control <- function(tail_completion = c("zero", "exp", "tau-zero"),
                             tail_tau = Inf,
                             maxit = 100,
                             epsilon = 1e-4,
                             pmin = 1e-5,
                             save_call = TRUE,
                             verbose = 0,
                             ...)
{
    dot_list <- list(...)
    ## for backward compatibility
    default_maxit <- missing(maxit)
    if ((isTRUE(dot_list$default_maxit) || default_maxit) &&
        !is.null(dot_list$max_it)) {
        maxit <- dot_list$max_it
    }
    default_epsilon <- missing(epsilon)
    if ((isTRUE(dot_list$default_epsilon) || default_epsilon) &&
        !is.null(dot_list$rel_tol)) {
        epsilon <- dot_list$rel_tol
    }
    if (is.integer(tail_completion)) {
        int_tail_completion <- tail_completion
    } else {
        tail_completion <- match.arg(tail_completion)
        tail_comps <- c("zero", "exp", "tau-zero")
        int_tail_completion <- match(tail_completion, tail_comps, nomatch = 1L)
    }
    ## return
    structure(
        list(
            tail_completion = int_tail_completion,
            tail_tau = {
                if (is.null(tail_tau)) {Inf} else {tail_tau}
            },
            maxit = maxit,
            epsilon = epsilon,
            pmin = pmin,
            save_call = save_call,
            verbose = verbose,
            default_maxit = default_maxit,
            default_epsilon = default_epsilon
        ),
    class = "cox_cure.control")
}


##' @rdname cox_cure
##'
##' @param start A numeric vector representing the initial values for the
##'     underlying model estimation procedure.  If \code{standardize} is
##'     \code{TRUE}, the specified initial values will be scaled internally to
##'     match the standardized data.  The default initial values depend on the
##'     specific models and based on the observed data.  If inappropriate
##'     initial values (in terms of length) are specified, the default values
##'     will be used.
##' @param offset A numeric vector specifying the offset term.  The length of
##'     the specified offset term should be equal to the sample size.
##' @param standardize A logical value specifying if each covariate should be
##'     standardized to have mean zero and standard deviation one internally for
##'     numerically stability and fair regularization.  The default value is
##'     \code{TRUE}.  The coefficient estimates will always be returned in
##'     original scales.
##' @export
cox_cure.mstep <- function(start = NULL,
                           offset = NULL,
                           maxit = 10,
                           epsilon = 1e-4,
                           standardize = TRUE,
                           ...)
{
    dot_list <- list(...)
    ## for backward compatibility
    default_maxit <- missing(maxit)
    if ((isTRUE(dot_list$default_maxit) || default_maxit) &&
        !is.null(dot_list$max_it)) {
        maxit <- dot_list$max_it
    }
    default_epsilon <- missing(epsilon)
    if ((isTRUE(dot_list$default_epsilon) || default_epsilon) &&
        !is.null(dot_list$rel_tol)) {
        epsilon <- dot_list$rel_tol
    }

    structure(
        list(
            start = null2num0(start),
            offset = null2num0(offset),
            maxit = maxit,
            epsilon = epsilon,
            standardize = standardize,
            default_maxit = default_maxit,
            default_epsilon = default_epsilon
        ),
        class = "cox_cure.mstep"
    )
}
