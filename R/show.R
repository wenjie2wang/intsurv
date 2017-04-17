## collation after class.R
##' @include class.R
NULL


##' Show Methods
##'
##' S4 class methods that display or summarize certain objects.
##'
##' \itemize{
##'
##' \item For \code{\link{intCox-class}} object, it prints out a brief summary
##'     of the fitted model.
##'
##' \item For \code{\link{intCox.summary-class}} object, it prints out summary
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
##' @aliases show,intCox-method
##' @importFrom stats printCoefmat
##' @export
setMethod(
    f = "show", signature = "intCox",
    definition = function(object)
    {
        ## function call
        cat("Call:\n", paste(deparse(Call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
        ## coefficient estimates
        cat("Coefficients of covariates: \n")
        printCoefmat(object@coefMat)
        ## return
        invisible(object)
    }
)


##' @rdname show-method
##' @aliases show,intCox.summary-method
##' @importFrom stats printCoefmat
##' @export
setMethod(
    f = "show", signature = "intCox.summary",
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
