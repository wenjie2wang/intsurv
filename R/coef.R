## collation after class.R
##' @include class.R
NULL


##' Estimated Coefficients of Covariates
##'
##' \code{coef,coxphx-method} is an S4 class method that extracts covariate
##' coefficient estimates from \code{\link{coxphx-class}} object from
##' function \code{\link{coxphx}}.
##'
##' @param object \code{\link{coxphx-class}} object.
##' @param ... Other arguments for future usage.
##' @return A named numeric vector.
##' @aliases coef,coxphx-method
##' @seealso
##' \code{\link{coxphx}} for fitting integrative Cox model;
##' \code{\link{summary,coxphx-method}} for summary of a fitted model.
##' @examples
##' ## See examples given in function coxphx.
##' @importFrom stats coef
##' @export
setMethod(
    f = "coef", signature = "coxphx",
    definition = function(object, ...)
    {
        object@estimates$beta[, "coef"]
    })
