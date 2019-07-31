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
##' model.
##'
##' @param object Object representing a fitted model.
##' @param ... Other arguments for future usage.
##'
##' @return A named numeric vector.
##' @importFrom stats coef
##' @method coef cox_cure
##' @export
coef.cox_cure <- function(object, ...)
{
    object$coef
}


##' Estimated Covariate Coefficients
##'
##' Extract the covariate coefficient estimates from a fitted Cox cure rate
##' model.
##'
##'
##' @param object Object representing a fitted model.
##' @param selection A character string for specifying the criterion for
##'     selection of coefficient estimates.  The available options are
##'     \code{"bic"} for selecting coefficient estimates based on BIC criterion
##'     and \code{"none"} for not performing selection and returning the whole
##'     solution path.
##' @param ... Other arguments for future usage.
##'
##' @return A named numeric vector.
##' @importFrom stats coef
##' @method coef cox_cure_net
##' @export
coef.cox_cure_net <- function(object, selection = c("bic", "none"), ...)
{
    object$coef[, which.min(object$model$bic)]
}
