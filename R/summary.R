## collation after class.R
##' @include class.R
NULL


##' Summary of a Fitted Model
##'
##' For \code{\link{intCox}} object, the function returns a
##' \code{\link{intCox.summary-class}} object whose slots include
##' \itemize{
##'     \item \code{call}: Function call of model fitting.
##'     \item \code{coefMat}: Estimated covariate coefficients.
##'     \item \code{logL}: Log-likelihood under observed data.
##' }
##'
##' @param object \code{\link{intCox-class}} object.
##' @param showCall A logic value with default \code{TRUE}, indicating whether
##'     \code{show} method prints out the original call information of
##'     \code{intCox}.  Set \code{FALSE} for a more concise printout.
##' @param ... Other arguments for future usage.
##'
##' @return \code{summary_intCox} class object.
##' @aliases summary,intCox-method
##' @examples
##' ## See examples given in function intCox
##' @seealso
##' \code{\link{intCox}} for model fitting;
##' \code{\link{coef,intCox}} for coefficient estimates.
##' @export
setMethod(
    f = "summary", signature = "intCox",
    definition = function(object, showCall = TRUE, ...)
    {
        Call <- object@call
        attr(Call, "show") <- showCall
        betaMat <- object@estimates$beta
        len_logL <- length(object@logL)
        new("intCox.summary",
            call = Call,
            coefMat = betaMat,
            logL = object@logL[len_logL])
    }
)
