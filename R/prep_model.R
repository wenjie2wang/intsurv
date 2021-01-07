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


### prepare design matrix and response from formula inputs for survival models
##' @importFrom stats model.matrix.default model.frame.default .getXlevels
prep_cure_model <- function(surv_formula, cure_formula,
                            obs_time, obs_event,
                            data, subset, contrasts = NULL,
                            eval_env = parent.frame())
{
    this_call <- call0 <- match.call(expand.dots = FALSE)
    if (missing(surv_formula)) {
        stop("The survival model formula cannot be missing.",
             call. = FALSE)
    }
    if (missing(cure_formula)) {
        stop("The cure model formula cannot be missing.",
             call. = FALSE)
    }
    if (missing(obs_time)) {
        stop("The survival times cannot be missing.",
             call. = FALSE)
    }
    if (missing(obs_time)) {
        stop("The event indicators cannot be missing.",
             call. = FALSE)
    }

    ## 1. process formula for survival model
    names(this_call)[names(this_call) == "surv_formula"] <- "formula"
    if (missing(data)) {
        this_call$data <- eval_env
    }
    matched_call <- match(c("formula", "data", "subset",
                            "obs_time", "obs_event"),
                          names(this_call), nomatch = 0L)
    this_call <- this_call[c(1L, matched_call)]
    ## drop unused levels in factors
    this_call$drop.unused.levels <- TRUE
    this_call$na.action <- na.pass
    this_call[[1L]] <- quote(stats::model.frame.default)
    mf <- eval(this_call, eval_env)
    surv_formula <- attr(mf, "terms")
    mf2 <- stats::model.frame.default(surv_formula, mf,
                                      na.action = na.fail)
    ## suppress warnings on not used contrasts
    suppressWarnings({
        mm <- stats::model.matrix(surv_formula, data = mf2,
                                  contrasts.arg = contrasts)
    })
    ## output: contrasts
    contrasts <- attr(mm, "contrasts")
    ## surv list
    surv_list <- list(
        obs_time = mf[["(obs_time)"]],
        obs_event = mf[["(obs_event)"]],
        xMat = mm,
        contrasts = contrasts,
        xlevels = stats::.getXlevels(attr(mf, "terms"), mf)
    )

    ## 2. process cure formula
    this_call <- call0
    names(this_call)[names(this_call) == "cure_formula"] <- "formula"
    if (missing(data)) {
        this_call$data <- eval_env
    }
    matched_call <- match(c("formula", "data", "subset"),
                          names(this_call), nomatch = 0L)
    this_call <- this_call[c(1L, matched_call)]
    ## drop unused levels in factors
    this_call$drop.unused.levels <- TRUE
    this_call$na.action <- na.fail
    this_call[[1L]] <- quote(stats::model.frame.default)
    mf <- eval(this_call, eval_env)
    cure_formula <- attr(mf, "terms")
    ## suppress warnings on not used contrasts
    suppressWarnings({
        mm <- stats::model.matrix.default(cure_formula, data = mf,
                                          contrasts.arg = contrasts)
    })
    ## output: contrasts
    contrasts <- attr(mm, "contrasts")
    ## cure list
    cure_list <- list(
        xMat = mm,
        contrasts = contrasts,
        xlevels = stats::.getXlevels(attr(mf, "terms"), mf)
    )

    ## return
    list(surv = surv_list,
         cure = cure_list)
}
