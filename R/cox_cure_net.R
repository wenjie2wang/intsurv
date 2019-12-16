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


##' Regularized Cox Cure Rate Model
##'
##' For right-censored data, fit a regularized Cox cure rate model through
##' elastic-net penalty following Masud et al. (2018), and Zou and Hastie
##' (2005).  For right-censored data with uncertain event status,
##' fit the regularized Cox cure model proposed by Wang et al. (2019+).  Without
##' regularization, the model reduces to the regular Cox cure rate model (Kuk
##' and Chen, 1992; Sy and Taylor, 2000)
##'
##' The model estimation procedure follows expectation maximization (EM)
##' algorithm.  Variable selection procedure through regularization by elastic
##' net penalty is developed based on cyclic coordinate descent and
##' majorization-minimization (MM) algorithm.
##'
##' @usage
##'
##' cox_cure_net(surv_formula, cure_formula, time, event,
##'              data, subset, contrasts = NULL,
##'              surv_lambda, surv_alpha = 1, surv_nlambda = 10,
##'              surv_lambda_min_ratio = 1e-1, surv_l1_penalty_factor,
##'              cure_lambda, cure_alpha = 1, cure_nlambda = 10,
##'              cure_lambda_min_ratio = 1e-1, cure_l1_penalty_factor,
##'              surv_start, cure_start,
##'              surv_standardize = TRUE, cure_standardize = TRUE,
##'              em_max_iter = 200, em_rel_tol = 1e-5,
##'              surv_max_iter = 10, surv_rel_tol = 1e-5,
##'              cure_max_iter = 10, cure_rel_tol = 1e-5,
##'              tail_completion = c("zero", "exp", "zero-tau"),
##'              tail_tau = NULL, pmin = 1e-5, early_stop = TRUE,
##'              verbose = FALSE, ...)
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
##' @param surv_lambda A numeric vector consists of non-negative values
##'     representing the tuning parameter sequence for the survival model part.
##' @param surv_alpha A number between 0 and 1 for tuning the elastic net
##'     penalty for the survival model part.  If it is one, the elastic penalty
##'     will reduce to the well-known lasso penalty.  If it is zero, the ridge
##'     penalty will be used.
##' @param surv_nlambda A positive number specifying the number of
##'     \code{surv_lambda} if \code{surv_lambda} is not specified.  The default
##'     value is 10.
##' @param surv_lambda_min_ratio The ratio of the minimum \code{surv_lambda} to
##'     the large enough \code{surv_lambda} that produces all-zero estimates on
##'     log scale.  The default value is \code{1e-1}.
##' @param surv_l1_penalty_factor A numeric vector that consists of positive
##'     numbers for penalty factors (or weights) on L1-norm for the coefficient
##'     estimate vector in the survival model part.  The penalty is applied to
##'     the coefficient estimate divided by the specified weights.  The
##'     specified weights are re-scaled internally so that their summation
##'     equals the length of coefficients.  If it is left unspecified, the
##'     weights are all set to be one.
##' @param cure_lambda A numeric vector consists of non-negative values
##'     representing the tuning parameter sequence for the cure model part.
##' @param cure_alpha A number between 0 and 1 for tuning the elastic net
##'     penalty for the cure model part.  If it is one, the elastic penalty will
##'     reduce to the well-known lasso penalty.  If it is zero, the ridge
##'     penalty will be used.
##' @param cure_nlambda A positive number specifying the number of
##'     \code{cure_lambda} if \code{cure_lambda} is not specified.  The default
##'     value is 10.
##' @param cure_lambda_min_ratio The ratio of the minimum \code{cure_lambda} to
##'     the large enough \code{cure_lambda} that produces all-zero estimates on
##'     log scale.  The default value is \code{1e-1}.
##' @param cure_l1_penalty_factor A numeric vector that consists of positive
##'     numbers for penalty factors (or weights) on L1-norm for the coefficient
##'     estimate vector in the cure model part.  The penalty is applied to the
##'     coefficient estimate divided by the specified weights.  The specified
##'     weights are re-scaled internally so that their summation equals the
##'     length of coefficients.  If it is left unspecified, the weights are all
##'     set to be one.
##' @param surv_start An optional numeric vector representing the starting
##'     values for the Cox model component.  If not specified, the starting
##'     values will be obtained from fitting a regular Cox model to events only.
##' @param cure_start An optional numeric vector representing the starting
##'     values for the logistic model component.  If not specified, the starting
##'     values will be obtained from fitting a regular logistic model to the
##'     non-missing event indicators.
##' @param surv_standardize A logical value specifying whether to standardize
##'     the covariates for the survival model part.  If \code{FALSE}, the
##'     covariates will be standardized internally to have mean zero and
##'     standard deviation one.
##' @param cure_standardize A logical value specifying whether to standardize
##'     the covariates for the cure rate model part.  If \code{TRUE} (by
##'     default), the covariates will be standardized internally to have mean
##'     zero and standard deviation one.
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
##'     The default value is \code{10} to encourage faster convergence.
##' @param surv_rel_tol A positive number specifying the tolerance that
##'     determines the convergence of the M-step related to the survival model
##'     component in terms of the convergence of the covariate coefficient
##'     estimates.  The tolerance is compared with the relative change between
##'     estimates from two consecutive iterations, which is measured by ratio of
##'     the L1-norm of their difference to the sum of their L1-norm.  The
##'     default value is \code{1e-5}.
##' @param cure_max_iter A positive integer specifying the maximum iteration
##'     number of the M-step routine related to the cure rate model component.
##'     The default value is \code{10} to encourage faster convergence.
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
##' \code{cox_cure_net} object for regular Cox cure rate model or
##' \code{cox_cure_net_uncer} object for Cox cure rate model with uncertain
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
##' Wang, W., Chen, K., Luo, C., & Yan, J. (2019+). Cox Cure Model with
##' Uncertain Event Status with application to a Suicide Risk
##' Study. \emph{Working in Progress}.
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
cox_cure_net <-
    function(surv_formula, cure_formula, time, event,
             data, subset, contrasts = NULL,
             surv_lambda,
             surv_alpha = 1,
             surv_nlambda = 10,
             surv_lambda_min_ratio = 1e-1,
             surv_l1_penalty_factor,
             cure_lambda,
             cure_alpha = 1,
             cure_nlambda = 10,
             cure_lambda_min_ratio = 1e-1,
             cure_l1_penalty_factor,
             surv_start,
             cure_start,
             surv_standardize = TRUE,
             cure_standardize = TRUE,
             em_max_iter = 200,
             em_rel_tol = 1e-5,
             surv_max_iter = 10,
             surv_rel_tol = 1e-5,
             cure_max_iter = 10,
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
    if (missing(surv_start)) {
        surv_start <- 0
    } else if (length(surv_start) != surv_x) {
        stop("The length of 'surv_start' is inappropriate.")
    }
    if (missing(cure_start)) {
        cure_start <- 0
    } else if (length(cure_start) != cure_x + as.integer(cure_intercept)) {
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

    ## lambda sequence
    if (missing(surv_lambda)) {
        surv_lambda <- 0
    } else {
        surv_nlambda <- 1
    }
    if (missing(cure_lambda)) {
        cure_lambda <- 0
    } else {
        cure_nlambda <- 1
    }

    ## penalty factor
    if (missing(surv_l1_penalty_factor)) {
        surv_l1_penalty_factor <- 0
    }
    if (missing(cure_l1_penalty_factor)) {
        cure_l1_penalty_factor <- 0
    }

    ## call the routine
    if (anyNA(obs_event)) {
        out <- coxph_cure_uncer_vs(
            time = obs_time,
            event = obs_event,
            cox_x = surv_x,
            cure_x = cure_x,
            cure_intercept = cure_intercept,
            cox_lambda = surv_lambda,
            cox_alpha = surv_alpha,
            cox_nlambda = surv_nlambda,
            cox_lambda_min_ratio = surv_lambda_min_ratio,
            cox_l1_penalty_factor = surv_l1_penalty_factor,
            cure_lambda = cure_lambda,
            cure_alpha = cure_alpha,
            cure_nlambda = cure_nlambda,
            cure_lambda_min_ratio = cure_lambda_min_ratio,
            cure_l1_penalty_factor = cure_l1_penalty_factor,
            cox_start = surv_start,
            cure_start = cure_start,
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
        class(out) <- "cox_cure_net_uncer"
    } else {
        out <- rcpp_coxph_cure_vs(
            time = obs_time,
            event = obs_event,
            cox_x = surv_x,
            cure_x = cure_x,
            cure_intercept = cure_intercept,
            cox_lambda = surv_lambda,
            cox_alpha = surv_alpha,
            cox_nlambda = surv_nlambda,
            cox_lambda_min_ratio = surv_lambda_min_ratio,
            cox_l1_penalty_factor = surv_l1_penalty_factor,
            cure_lambda = cure_lambda,
            cure_alpha = cure_alpha,
            cure_nlambda = cure_nlambda,
            cure_lambda_min_ratio = cure_lambda_min_ratio,
            cure_l1_penalty_factor = cure_l1_penalty_factor,
            cox_start = surv_start,
            cure_start = cure_start,
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
                  paste0("z", seq_len(ncol(out$cure_coef[- 1L]))))
            } else {
                paste0("z", seq_len(ncol(out$cure_coef)))
            }
    }
    colnames(out$surv_en_coef) <- colnames(out$surv_coef)
    colnames(out$cure_en_coef) <- colnames(out$cure_coef)
    ## add function call
    out$call <- call0
    ## return
    out
}


##' @rdname cox_cure_net
##'
##' @usage
##'
##' cox_cure_net.fit(surv_x, cure_x, time, event, cure_intercept = TRUE,
##'                  surv_lambda, surv_alpha = 1, surv_nlambda = 10,
##'                  surv_lambda_min_ratio = 1e-1, surv_l1_penalty_factor,
##'                  cure_lambda, cure_alpha = 1, cure_nlambda = 10,
##'                  cure_lambda_min_ratio = 1e-1, cure_l1_penalty_factor,
##'                  surv_start, cure_start,
##'                  surv_standardize = TRUE, cure_standardize = TRUE,
##'                  em_max_iter = 200, em_rel_tol = 1e-5,
##'                  surv_max_iter = 10, surv_rel_tol = 1e-5,
##'                  cure_max_iter = 10, cure_rel_tol = 1e-5,
##'                  tail_completion = c("zero", "exp", "zero-tau"),
##'                  tail_tau = NULL, pmin = 1e-5, early_stop = TRUE,
##'                  verbose = FALSE, ...)
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
cox_cure_net.fit <-
    function(surv_x, cure_x, time, event,
             cure_intercept = TRUE,
             surv_lambda,
             surv_alpha = 1,
             surv_nlambda = 10,
             surv_lambda_min_ratio = 1e-1,
             surv_l1_penalty_factor,
             cure_lambda,
             cure_alpha = 1,
             cure_nlambda = 10,
             cure_lambda_min_ratio = 1e-1,
             cure_l1_penalty_factor,
             surv_start,
             cure_start,
             surv_standardize = TRUE,
             cure_standardize = TRUE,
             em_max_iter = 200,
             em_rel_tol = 1e-5,
             surv_max_iter = 10,
             surv_rel_tol = 1e-5,
             cure_max_iter = 10,
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
    if (missing(surv_start)) {
        surv_start <- 0
    } else if (length(surv_start) != surv_x) {
        stop("The length of 'surv_start' is inappropriate.")
    }
    if (missing(cure_start)) {
        cure_start <- 0
    } else if (length(cure_start) != cure_x + as.integer(cure_intercept)) {
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

    ## lambda sequence
    if (missing(surv_lambda)) {
        surv_lambda <- 0
    } else {
        surv_nlambda <- 1
    }
    if (missing(cure_lambda)) {
        cure_lambda <- 0
    } else {
        cure_nlambda <- 1
    }

    ## penalty factor
    if (missing(surv_l1_penalty_factor)) {
        surv_l1_penalty_factor <- 0
    }
    if (missing(cure_l1_penalty_factor)) {
        cure_l1_penalty_factor <- 0
    }

    ## call the routine
    if (anyNA(event)) {
        out <- coxph_cure_uncer_vs(
            time = time,
            event = event,
            cox_x = surv_x,
            cure_x = cure_x,
            cure_intercept = cure_intercept,
            cox_lambda = surv_lambda,
            cox_alpha = surv_alpha,
            cox_nlambda = surv_nlambda,
            cox_lambda_min_ratio = surv_lambda_min_ratio,
            cox_l1_penalty_factor = surv_l1_penalty_factor,
            cure_lambda = cure_lambda,
            cure_alpha = cure_alpha,
            cure_nlambda = cure_nlambda,
            cure_lambda_min_ratio = cure_lambda_min_ratio,
            cure_l1_penalty_factor = cure_l1_penalty_factor,
            cox_start = surv_start,
            cure_start = cure_start,
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
        class(out) <- "cox_cure_net_uncer"
    } else {
        out <- rcpp_coxph_cure_vs(
            time = time,
            event = event,
            cox_x = surv_x,
            cure_x = cure_x,
            cure_intercept = cure_intercept,
            cox_lambda = surv_lambda,
            cox_alpha = surv_alpha,
            cox_nlambda = surv_nlambda,
            cox_lambda_min_ratio = surv_lambda_min_ratio,
            cox_l1_penalty_factor = surv_l1_penalty_factor,
            cure_lambda = cure_lambda,
            cure_alpha = cure_alpha,
            cure_nlambda = cure_nlambda,
            cure_lambda_min_ratio = cure_lambda_min_ratio,
            cure_l1_penalty_factor = cure_l1_penalty_factor,
            cox_start = surv_start,
            cure_start = cure_start,
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
                  paste0("z", seq_len(ncol(out$cure_coef[- 1L]))))
            } else {
                paste0("z", seq_len(ncol(out$cure_coef)))
            }
    }
    colnames(out$surv_en_coef) <- colnames(out$surv_coef)
    colnames(out$cure_en_coef) <- colnames(out$cure_coef)

    ## add function call
    out$call <- this_call
    ## return
    out
}
