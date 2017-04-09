## collation after class.R
##' @include class.R
NULL


##' FIXME
setMethod(
    f = "summary", signature = "intCox",
    definition = function(object, showCall = TRUE, ...)
    {
        Call <- object@call
        attr(Call, "show") <- showCall
        betaMat <- object@estimates$beta
        new("summary_intCox",
            call = Call,
            coefMat = betaMat)
    }
)
