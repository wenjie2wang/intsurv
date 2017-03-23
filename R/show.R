## collation after class.R
##' @include class.R
NULL


##' Show an object.
##'
##' An S4 class generic function that displays certain object.
##'
##' \itemize{
##'   \item For \code{\link{summaryCoxEm-class}} object,
##'       it prints out summary of the fitted model.
##' }
##'
##' @param object An object used to dispatch a method.
##' @name show-method
##' @importFrom methods show
NULL


##' @rdname show-method
##' @aliases show,summaryCoxEm-method
##' @importFrom stats printCoefmat
##' @export
setMethod(
    f = "show", signature = "summaryCoxEm",
    definition = function(object) {
        ## if (attr(object@call, "show")) {
        ##     Call <- object@call
        ##     attr(Call, "show") <- NULL
        ##     cat("Call: \n")
        ##     print(Call)
        ##     cat("\n")
        ## }
        cat("Coefficients of covariates: \n")
        printCoefmat(object@coefMat)
    }
)
