## collation after class.R
##' @include class.R
NULL


##' Summary of a Fitted Model
##'
##' For \code{\link{coxphx}} object, the function returns a
##' \code{\link{coxphx.summary-class}} object whose slots include
##' \itemize{
##'     \item \code{call}: Function call of model fitting.
##'     \item \code{coefMat}: Estimated covariate coefficients.
##'     \item \code{logL}: Log-likelihood under observed data.
##' }
##'
##' @param object \code{\link{coxphx-class}} object.
##' @param showCall A logic value with default \code{TRUE}, indicating whether
##'     \code{show} method prints out the original call information of
##'     \code{coxphx}.  Set \code{FALSE} for a more concise printout.
##' @param ... Other arguments for future usage.
##'
##' @return \code{coxphx.summary-class} class object.
##' @aliases summary,coxphx-method
##' @examples
##' ## See examples given in function coxphx
##' @seealso
##' \code{\link{coxphx}} for model fitting;
##' \code{\link{coef,coxphx}} for coefficient estimates.
##' @export
setMethod(
    f = "summary", signature = "coxphx",
    definition = function(object, showCall = TRUE, ...)
    {
        Call <- object@call
        attr(Call, "show") <- showCall
        betaMat <- object@estimates$beta
        len_logL <- length(object@logL)
        new("coxphx.summary",
            call = Call,
            coefMat = betaMat,
            logL = object@logL[len_logL])
    }
)
