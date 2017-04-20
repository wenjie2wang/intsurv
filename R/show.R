################################################################################
##
##   R package intsurv by Wenjie Wang, Kun Chen and Jun Yan
##   Copyright (C) 2017
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


##' Show Methods
##'
##' S4 class methods that display or summarize certain objects.
##'
##' \itemize{
##'
##' \item For \code{\link{coxphx-class}} object, it prints out a brief summary
##'     of the fitted model.
##'
##' \item For \code{\link{coxphx.summary-class}} object, it prints out summary
##'       of a fitted model.
##'
##' }
##'
##' @param object An object used to dispatch a method.
##' @return object input (invisibly).
##' @name show-method
##' @importFrom methods show
NULL


##' @rdname show-method
##' @aliases show,coxphx-method
##' @importFrom stats printCoefmat
##' @export
setMethod(
    f = "show", signature = "coxphx",
    definition = function(object)
    {
        ## function call
        cat("Call:\n", paste(deparse(object@call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
        ## coefficient estimates
        cat("Coefficients of covariates: \n")
        printCoefmat(object@estimates$beta)
        ## return
        invisible(object)
    }
)


##' @rdname show-method
##' @aliases show,coxphx.summary-method
##' @importFrom stats printCoefmat
##' @export
setMethod(
    f = "show", signature = "coxphx.summary",
    definition = function(object)
    {
        ## function call
        if (attr(object@call, "show")) {
            Call <- object@call
            attr(Call, "show") <- NULL
            cat("Call:\n",
                paste(deparse(Call), sep = "\n",
                      collapse = "\n"),
                "\n\n", sep = "")
        }
        ## coefficient estimates
        cat("Coefficients of covariates: \n")
        printCoefmat(object@coefMat)
        ## log-likelihood after convergence
        cat("\nLog-likelihood under observed data:", object@logL, "\n")
        ## return
        invisible(object)
    }
)
