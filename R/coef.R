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


## collation after class.R
##' @include class.R
NULL


##' Estimated Covariates Coefficients
##'
##' \code{coef,iCoxph-method} is an S4 class method that extracts covariate
##' coefficient estimates from \code{\link{iCoxph-class}} object from function
##' \code{\link{iCoxph}}.
##'
##' @param object \code{\link{iCoxph-class}} object.
##' @param ... Other arguments for future usage.
##' @return A named numeric vector.
##' @aliases coef,iCoxph-method coef.iCoxph
##' @seealso
##' \code{\link{iCoxph}} for fitting integrative Cox model;
##' \code{\link{summary,iCoxph-method}} for summary of a fitted model.
##' @examples
##' ## See examples of function iCoxph.
##' @importFrom stats coef
##' @export
setMethod(
    f = "coef", signature = "iCoxph",
    definition = function(object, ...)
    {
        object@estimates$beta[, "coef"]
    })


##' Estimated Covariate Coefficients
##'
##' Extract the covariate coefficient estimates from a fitted Cox cure rate
##' model with possible uncertain event status.
##'
##' @param object Object representing a fitted model.
##' @param part A character string specifying the coefficient estimates from a
##'     particular model part.  The available options are \code{"both"} for the
##'     coefficient estimates from both model parts, \code{"survival"} for the
##'     coefficient estimates from the survival model part, and \code{"cure"}
##'     for the coefficient estimates from the cure rate model part.  The
##'     default option is \code{"all"}.
##' @param ... Other arguments for future usage.
##'
##' @return If \code{part = "both"}, this function returns a list that consists
##'     of the following named elements \itemize{
##'
##' \item \code{surv}: the coefficient estimates of survival model part.
##'
##' \item \code{cure}: the coefficient estimates of cure rate model part.
##'
##' }
##'
##' Otherwise, a named numeric vector representing the coefficient estimates of
##' the specified model part will be returned.
##'
##' @importFrom stats coef
##' @method coef cox_cure
##' @export
coef.cox_cure <- function(object, part = c("both", "survival", "cure"), ...)
{
    part <- match.arg(part)
    if (part != "both") {
        coef_name <- sprintf(
            "%s_coef", switch(part, "survival" = "surv", "cure" = "cure")
        )
        return(object[[coef_name]])
    }
    ## else return both
    list(surv = object$surv_coef,
         cure = object$cure_coef)
}


##' @rdname coef.cox_cure
##' @method coef cox_cure_mar
##' @export
coef.cox_cure_mar <- coef.cox_cure


##' Estimated Covariate Coefficients
##'
##' Extract the covariate coefficient estimates from a solution path of
##' regularized Cox cure rate model.
##'
##' @param object Object representing a fitted solution path.
##' @param selection A character string for specifying the criterion for
##'     selection of coefficient estimates.  The available options are
##'     \code{"bic1"} for selecting coefficient estimates by regular BIC
##'     criterion based on the number of observations, \code{"bic2"} for a
##'     variant BIC criterion based on the effective sample size, and
##'     \code{"all"} for returning the whole solution path.  See
##'     \code{\link{BIC.cox_cure_net}} for details of selection by BIC.
##' @param ... Other arguments for future usage.
##'
##' @return A list that consists of the following named elements:
##' \itemize{
##'
##' \item \code{surv}: the selected coefficient estimates of survival model
##'     part.
##'
##' \item \code{cure}: the selected coefficient estimates of cure rate model
##'     part.
##'
##' }
##'
##' @examples
##' ## see examples of function `cox_cure_net`
##'
##' @importFrom stats coef
##' @method coef cox_cure_net
##' @export
coef.cox_cure_net <- function(object,
                              selection = c("bic1", "bic2", "all"), ...)
{
    bic <- match.arg(selection)
    if (bic != "all") {
        idx <- which.min(object$model[[bic]])
        out <- list(surv = object$surv_coef[idx, ],
                    cure = object$cure_coef[idx, ])
    } else {
        out <- list(surv = object$surv_coef,
                    cure = object$cure_coef)
    }
    ## return
    out
}


##' @rdname coef.cox_cure_net
##' @importFrom stats coef
##' @method coef cox_cure_net_mar
##' @export
coef.cox_cure_net_mar <- coef.cox_cure_net
