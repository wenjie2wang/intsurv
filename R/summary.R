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


##' Summary of a Fitted Model
##'
##' For \code{\link{iCoxph}} object, the function returns a
##' \code{\link{iCoxph.summary-class}} object whose slots include
##' \itemize{
##'     \item \code{call}: Function call of model fitting.
##'     \item \code{coefMat}: Estimated covariate coefficients.
##'     \item \code{logL}: Log-likelihood under observed data.
##' }
##'
##' @param object \code{\link{iCoxph-class}} object.
##' @param showCall A logic value with default \code{TRUE}, indicating whether
##'     \code{show} method prints out the original call information of
##'     \code{iCoxph}.  Set \code{FALSE} for a more concise printout.
##' @param ... Other arguments for future usage.
##'
##' @return \code{iCoxph.summary-class} class object.
##' @aliases summary,iCoxph-method
##' @examples
##' ## See examples given in function iCoxph
##' @seealso
##' \code{\link{iCoxph}} for model fitting;
##' \code{\link{coef,iCoxph-method}} for coefficient estimates.
##' @export
setMethod(
    f = "summary", signature = "iCoxph",
    definition = function(object, showCall = TRUE, ...)
    {
        Call <- object@call
        attr(Call, "show") <- showCall
        betaMat <- object@estimates$beta
        len_logL <- length(object@logL)
        new("iCoxph.summary",
            call = Call,
            coefMat = betaMat,
            logL = object@logL[len_logL])
    }
)
