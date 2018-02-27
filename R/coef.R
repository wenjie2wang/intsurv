################################################################################
##
##   R package intsurv by Wenjie Wang, Kun Chen and Jun Yan
##   Copyright (C) 2017-2018
##
##   This file is part of the R package intsurv.
##
##   The R package intsurv is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package intsurv is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################


## collation after class.R
##' @include class.R
NULL


##' Estimated Coefficients of Covariates
##'
##' \code{coef,iCoxph-method} is an S4 class method that extracts covariate
##' coefficient estimates from \code{\link{iCoxph-class}} object from
##' function \code{\link{iCoxph}}.
##'
##' @param object \code{\link{iCoxph-class}} object.
##' @param ... Other arguments for future usage.
##' @return A named numeric vector.
##' @aliases coef,iCoxph-method
##' @seealso
##' \code{\link{iCoxph}} for fitting integrative Cox model;
##' \code{\link{summary,iCoxph-method}} for summary of a fitted model.
##' @examples
##' ## See examples given in function iCoxph.
##' @importFrom stats coef
##' @export
setMethod(
    f = "coef", signature = "iCoxph",
    definition = function(object, ...)
    {
        object@estimates$beta[, "coef"]
    })
