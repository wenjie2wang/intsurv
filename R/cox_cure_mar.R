##' Cox Cure Rate Model For Uncertain Event Indicators
##'
##' For right-censored data with uncertain event indicators, fit the Cox cure
##' rate model proposed by Wang et al. (2020).
##'
##' @inheritParams cox_cure
##'
##' @param mar_formula A formula specifying a logistic model for modeling the
##'     missingness of the event indicators under the missing-at-random (MAR)
##'     assumption.
##' @param mar_control A list of control parameters for the MAR model part.
##' @param event A numeric vector representing the event indicators, where
##'     \code{NA}'s are allowed and represent uncertain event status.
##'
##' @return A \code{cox_cure_mar} object for the fitted Cox cure rate model with
##'     uncertain event indicators.
##'
##' @references
##'
##' Wang, W., Luo, C., Aseltine, R. H., Wang, F., Yan, J., & Chen, K. (2020).
##' Suicide Risk Modeling with Uncertain Diagnostic Records. \emph{arXiv
##' preprint arXiv:2009.02597}.
##'
##' @export
cox_cure_mar <- function(surv_formula,
                         cure_formula,
                         mar_formula,
                         time,
                         event,
                         data,
                         subset,
                         contrasts = NULL,
                         bootstrap = 0,
                         control = cox_cure.control(),
                         surv_control = intsurv.control(),
                         cure_control = intsurv.control(),
                         mar_control = intsurv.control(),
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
    if (! inherits(mar_control, "intsurv.control")) {
        mar_control <- do.call(intsurv.control, mar_control)
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
    names(mar_control) <- paste0("mar_", names(mar_control))
    call_list <- c(control, surv_control, cure_control, mar_control)
    ## get design matrix and responses
    call_list$time <- model_list$surv$time
    call_list$event <- model_list$surv$event
    call_list$surv_x <- model_list$surv$x
    call_list$cure_x <- model_list$cure$x
    call_list$mar_x <- model_list$mar$x
    if (! is.null(model_list$surv$offset)) {
        call_list$surv_offset <- model_list$surv$offset
    }
    if (! is.null(model_list$cure$offset)) {
        call_list$cure_offset <- model_list$cure$offset
    }
    if (! is.null(model_list$mar$offset)) {
        call_list$mar_offset <- model_list$mar$offset
    }
    call_list$bootstrap <- as.integer(bootstrap)
    ## for the surv part, remove the possible intercept term
    surv_is_intercept <- colnames(call_list$surv_x) == "(Intercept)"
    surv_has_intercept <- any(surv_is_intercept)
    if ((ncol(call_list$surv_x) - as.integer(surv_has_intercept)) == 0L) {
        stop("No covariate is specified in 'surv_formula'.")
    }
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
    mar_is_intercept <- colnames(call_list$mar_x) == "(Intercept)"
    mar_has_intercept <- any(mar_is_intercept)
    mar_only_intercept <- all(mar_is_intercept)
    call_list$mar_standardize <- ! mar_only_intercept
    if (mar_has_intercept) {
        call_list$mar_x <- call_list$mar_x[, which(! mar_is_intercept),
                                           drop = FALSE]
        call_list$mar_intercept <- TRUE
    } else {
        call_list$mar_intercept <- FALSE
    }
    ## check event
    if (all(call_list$event[! is.na(call_list$event)] < 1)) {
        stop("No event can be found.")
    }
    ## start values
    call_list$surv_start <- null2num0(call_list$surv_start)
    call_list$cure_start <- null2num0(call_list$cure_start)
    call_list$mar_start <- null2num0(call_list$mar_start)
    ## offset terms
    call_list$surv_offset <- null2num0(call_list$surv_offset)
    call_list$cure_offset <- null2num0(call_list$cure_offset)
    call_list$mar_offset <- null2num0(call_list$mar_offset)
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
    is_arg_valid <- names(call_list) %in%
        names(formals(rcpp_coxph_cure_mar))
    call_list <- call_list[is_arg_valid]
    out <- do.call(rcpp_coxph_cure_mar, call_list)
    ## add class
    class(out) <- "cox_cure_mar"
    ## add bootstrap se if available
    if (bootstrap <= 0) {
        out$bootstrap$surv_coef_mat <-
            out$bootstrap$cure_coef_mat <-
                out$bootstrap$mar_coef_mat <- NULL
    } else {
        out$bootstrap$surv_coef_se <- apply(out$bootstrap$surv_coef_mat,
                                            2L, se_interQ)
        out$bootstrap$cure_coef_se <- apply(out$bootstrap$cure_coef_mat,
                                            2L, se_interQ)
        out$bootstrap$mar_coef_se <- apply(out$bootstrap$mar_coef_mat,
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
    if (! is.null(mar_var_names <- colnames(call_list$mar_x))) {
        names(out$mar_coef) <- c(
        { if (call_list$mar_intercept) "(Intercept)" else NULL },
        mar_var_names
        )
    } else {
        names(out$mar_coef) <-
            if (call_list$mar_intercept) {
                c("(Intercept)",
                  paste0("v", seq_along(out$mar_coef[- 1L])))
            } else {
                paste0("v", seq_along(out$mar_coef))
            }
    }
    if (control$save_call) {
        ## add function call
        out$call <- call0
    }
    ## return
    out
}


##' @rdname cox_cure_mar
##'
##' @param mar_x A numeric matrix representing the design matrices for the MAR
##'     model part.
##' @param mar_intercept A logical value specifying whether to add an intercept
##'     term to the MAR model part.  If \code{TRUE} by default, an intercept
##'     term is included.
##'
##' @export
cox_cure_mar.fit <- function(surv_x,
                             cure_x,
                             mar_x,
                             time,
                             event,
                             cure_intercept = TRUE,
                             mar_intercept = TRUE,
                             bootstrap = 0,
                             control = cox_cure.control(),
                             surv_control = intsurv.control(),
                             cure_control = intsurv.control(),
                             mar_control = intsurv.control(),
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
    if (! inherits(mar_control, "intsurv.control")) {
        mar_control <- do.call(intsurv.control, mar_control)
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
    names(mar_control) <- paste0("mar_", names(mar_control))
    call_list <- c(control, surv_control, cure_control, mar_control)
    call_list$time <- time
    call_list$event <- event
    call_list$surv_x <- surv_x
    call_list$cure_x <- cure_x
    call_list$mar_x <- mar_x
    call_list$bootstrap <- bootstrap
    call_list$cure_intercept <- cure_intercept
    call_list$mar_intercept <- mar_intercept
    ## start values
    call_list$surv_start <- null2num0(call_list$surv_start)
    call_list$cure_start <- null2num0(call_list$cure_start)
    call_list$mar_start <- null2num0(call_list$mar_start)
    ## offset terms
    call_list$surv_offset <- null2num0(call_list$surv_offset)
    call_list$cure_offset <- null2num0(call_list$cure_offset)
    call_list$mar_offset <- null2num0(call_list$mar_offset)
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
    is_arg_valid <- names(call_list) %in%
        names(formals(rcpp_coxph_cure_mar))
    call_list <- call_list[is_arg_valid]
    out <- do.call(rcpp_coxph_cure_mar, call_list)
    ## add class
    class(out) <- "cox_cure_mar"
    ## add bootstrap se if available
    if (bootstrap <= 0) {
        out$bootstrap$surv_coef_mat <-
            out$bootstrap$cure_coef_mat <-
                out$bootstrap$mar_coef_mat <- NULL
    } else {
        out$bootstrap$surv_coef_se <- apply(out$bootstrap$surv_coef_mat,
                                            2L, se_interQ)
        out$bootstrap$cure_coef_se <- apply(out$bootstrap$cure_coef_mat,
                                            2L, se_interQ)
        out$bootstrap$mar_coef_se <- apply(out$bootstrap$mar_coef_mat,
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
    if (! is.null(mar_var_names <- colnames(mar_x))) {
        names(out$mar_coef) <- c({if (mar_intercept) "(Intercept)" else NULL},
                                 mar_var_names)
    } else {
        names(out$mar_coef) <-
            if (mar_intercept) {
                c("(Intercept)",
                  paste0("v", seq_along(out$mar_coef[- 1L])))
            } else {
                paste0("v", seq_along(out$mar_coef))
            }
    }
    if (control$save_call) {
        ## add function call
        out$call <- call0
    }
    ## return
    out
}


##' @rdname cox_cure_mar
##' @export
cox_cure_net_mar <- function(surv_formula,
                             cure_formula,
                             mar_formula,
                             time,
                             event,
                             data,
                             subset,
                             contrasts = NULL,
                             cv_nfolds = 0,
                             control = cox_cure.control(),
                             surv_control = cox_cure_net.control(),
                             cure_control = cox_cure_net.control(),
                             mar_control = cox_cure_net.control(),
                             ...)
{
    ## controls
    if (! inherits(control, "cox_cure.control")) {
        control <- do.call(cox_cure.control, control)
    }
    if (! inherits(surv_control, "cox_cure_net.control")) {
        surv_control <- do.call(cox_cure_net.control, surv_control)
    }
    if (! inherits(cure_control, "cox_cure_net.control")) {
        cure_control <- do.call(cox_cure_net.control, cure_control)
    }
    if (! inherits(mar_control, "cox_cure_net.control")) {
        mar_control <- do.call(cox_mar_net.control, mar_control)
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
    names(mar_control) <- paste0("mar_", names(mar_control))
    call_list <- c(control, surv_control, cure_control, mar_control)
    ## get design matrix and responses
    call_list$time <- model_list$surv$time
    call_list$event <- model_list$surv$event
    call_list$surv_x <- model_list$surv$x
    call_list$cure_x <- model_list$cure$x
    call_list$mar_x <- model_list$mar$x
    if (! is.null(model_list$surv$offset)) {
        call_list$surv_offset <- model_list$surv$offset
    }
    if (! is.null(model_list$cure$offset)) {
        call_list$cure_offset <- model_list$cure$offset
    }
    if (! is.null(model_list$mar$offset)) {
        call_list$mar_offset <- model_list$mar$offset
    }
    call_list$cv_nfolds <- as.integer(cv_nfolds)
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
    cure_standardize <- ! cure_only_intercept
    ## remove the possible intercept term
    if (cure_has_intercept) {
        call_list$cure_x <- call_list$cure_x[, which(! cure_is_intercept),
                                             drop = FALSE]
        call_list$cure_intercept <- TRUE
    } else {
        call_list$cure_intercept <- FALSE
    }
    mar_is_intercept <- colnames(call_list$mar_x) == "(Intercept)"
    mar_has_intercept <- any(mar_is_intercept)
    mar_only_intercept <- all(mar_is_intercept)
    mar_standardize <- ! mar_only_intercept
    ## remove the possible intercept term
    if (mar_has_intercept) {
        call_list$mar_x <- call_list$mar_x[, which(! mar_is_intercept),
                                           drop = FALSE]
        call_list$mar_intercept <- TRUE
    } else {
        call_list$mar_intercept <- FALSE
    }
    ## check event
    if (all(call_list$event[! is.na(call_list$event)] < 1)) {
        stop("No event can be found.")
    }
    ## start values
    call_list$surv_start <- null2num0(call_list$surv_start)
    call_list$cure_start <- null2num0(call_list$cure_start)
    call_list$mar_start <- null2num0(call_list$mar_start)
    ## offset terms
    call_list$surv_offset <- null2num0(call_list$surv_offset)
    call_list$cure_offset <- null2num0(call_list$cure_offset)
    call_list$mar_offset <- null2num0(call_list$mar_offset)
    ## more checks if tail completion after a specified tau
    if (call_list$tail_completion == 3L) {
        is_tau_small <- with(call_list,
                             tail_tau < max(time[! is.na(event) & event > 0]))
        if (is_tau_small) {
            stop("The specified 'tail_tau' cannot be less than",
                 " the largest event time.")
        }
    }
    is_arg_valid <- names(call_list) %in%
        names(formals(rcpp_coxph_cure_mar_vs))
    call_list <- call_list[is_arg_valid]
    out <- do.call(rcpp_coxph_cure_mar_vs, call_list)
    ## add class
    class(out) <- "cox_cure_net_mar"
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
    if (! is.null(mar_var_names <- colnames(call_list$mar_x))) {
        colnames(out$mar_coef) <- c(
        { if (call_list$mar_intercept) "(Intercept)" else NULL },
        mar_var_names
        )
    } else {
        colnames(out$mar_coef) <-
            if (call_list$mar_intercept) {
                c("(Intercept)",
                  paste0("v", seq_along(out$mar_coef[- 1L])))
            } else {
                paste0("v", seq_along(out$mar_coef))
            }
    }
    if (control$save_call) {
        ## add function call
        out$call <- call0
    }
    ## return
    out
}


##' @rdname cox_cure_mar
##' @export
cox_cure_net_mar.fit <- function(surv_x,
                                 cure_x,
                                 mar_x,
                                 time,
                                 event,
                                 cure_intercept = TRUE,
                                 mar_intercept = TRUE,
                                 cv_nfolds = 0,
                                 control = cox_cure.control(),
                                 surv_control = cox_cure_net.control(),
                                 cure_control = cox_cure_net.control(),
                                 mar_control = cox_cure_net.control(),
                                 ...)
{
    ## controls
    if (! inherits(control, "cox_cure.control")) {
        control <- do.call(cox_cure.control, control)
    }
    if (! inherits(surv_control, "cox_cure_net.control")) {
        surv_control <- do.call(cox_cure_net.control, surv_control)
    }
    if (! inherits(cure_control, "cox_cure_net.control")) {
        cure_control <- do.call(cox_cure_net.control, cure_control)
    }
    if (! inherits(mar_control, "cox_cure_net.control")) {
        mar_control <- do.call(cox_mar_net.control, mar_control)
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
    names(mar_control) <- paste0("mar_", names(mar_control))
    call_list <- c(control, surv_control, cure_control, mar_control)
    call_list$time <- time
    call_list$event <- event
    call_list$surv_x <- surv_x
    call_list$cure_x <- cure_x
    call_list$cv_nfolds <- cv_nfolds
    call_list$cure_intercept <- cure_intercept
    call_list$mar_intercept <- mar_intercept
    ## start values
    call_list$surv_start <- null2num0(call_list$surv_start)
    call_list$cure_start <- null2num0(call_list$cure_start)
    call_list$mar_start <- null2num0(call_list$mar_start)
    ## offset terms
    call_list$surv_offset <- null2num0(call_list$surv_offset)
    call_list$cure_offset <- null2num0(call_list$cure_offset)
    call_list$mar_offset <- null2num0(call_list$mar_offset)
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
    is_arg_valid <- names(call_list) %in%
        names(formals(rcpp_coxph_cure_mar_vs))
    call_list <- call_list[is_arg_valid]
    out <- do.call(rcpp_coxph_cure_mar_vs, call_list)
    ## add class
    class(out) <- "cox_cure_net_mar"
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
    if (! is.null(mar_var_names <- colnames(mar_x))) {
        colnames(out$mar_coef) <- c(
        {if (mar_intercept) "(Intercept)" else NULL},
        colnames(mar_x)
        )
    } else {
        colnames(out$mar_coef) <-
            if (mar_intercept) {
                c("(Intercept)",
                  paste0("v", seq_len(ncol(out$mar_coef) - 1L)))
            } else {
                paste0("v", seq_len(ncol(out$mar_coef)))
            }
    }
    if (control$save_call) {
        ## add function call
        out$call <- call0
    }
    ## return
    out
}
