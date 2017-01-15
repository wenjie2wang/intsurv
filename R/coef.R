## collation after class.R
#' @include class.R
NULL


#' Estimated Coefficients of Covariates
#'
#' \code{coef,coxEm-method} is a S4 class method that extracts
#' estimated coefficients of covariates from
#' \code{\link{coxEm-class}} object produced by
#' function \code{\link{coxEm}}.
#'
#' @param object \code{\link{coxEm-class}} object.
#' @param ... Other arguments for future usage.
#' @return A named numeric vector.
#' @aliases coef,coxEm-method
#' @seealso
#' \code{\link{coxEm}} for model fitting;
#' \code{\link{confint,coxEm-method}} for confidence intervals
#' for covariate coefficients;
#' \code{\link{summary,coxEm-method}} for summary of a fitted model.
#' @examples
#' ## See examples given in function coxEm.
#' @importFrom stats coef
#' @export
setMethod(f = "coef", signature = "coxEm",
          definition = function(object, ...) {
              object@estimates$beta
          })


