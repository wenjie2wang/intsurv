## collation after class.R
##' @include class.R
NULL


##' FIXME
setMethod(
    f = "summary", signature = "coxEm",
    definition = function(object, showCall = TRUE, ...)
    {
        Call <- object@call
        attr(Call, "show") <- showCall
        betaMat <- object@estimates$beta
        new("summaryCoxEm",
            call = Call,
            coefMat = betaMat)
    }
)
