## collation after class.R
##' @include class.R
NULL


##' Estimated Coefficients of Covariates
##'
##' \code{coef,intCox-method} is an S4 class method that extracts covariate
##' coefficient estimates from \code{\link{intCox-class}} object from
##' function \code{\link{intCox}}.
##'
##' @param object \code{\link{intCox-class}} object.
##' @param ... Other arguments for future usage.
##' @return A named numeric vector.
##' @aliases coef,intCox-method
##' @seealso
##' \code{\link{intCox}} for fitting integrative Cox model;
##' \code{\link{summary,intCox-method}} for summary of a fitted model.
##' @examples
##' ## See examples given in function intCox.
##' @importFrom stats coef
##' @export
setMethod(
    f = "coef", signature = "intCox",
    definition = function(object, ...)
    {
        object@estimates$beta[, "coef"]
    })
