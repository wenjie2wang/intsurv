#' Formula Response for survival data with uncertain records
#'
#' \code{Surve} is an S3 class that represents
#' formula response for survival data with
#' uncertain records.
#' The last letter 'e' in 'Surve' represents 'EM'.
#'
#' @param ID Identificator of each subject.
#' @param time Time of reccurence event or censoring.
#' @param event The status indicator, 0 = censored, 1 = event.
#' @aliases Surve
#' @export
Surve <- function (ID, time, event) {
    ## inpDat <- data.frame(ID, time, event)
    ## FIXME: function to check the date structure of
    ## uncertained records
    ## dat <- checkSurve(inpDat)
    ## outDat <- with(dat, as.matrix(cbind(ID, time, event)))
    outDat <- as.matrix(cbind(ID, time, event))
    ## attr(outDat, "ID") <- attr(dat, "ID")
    oldClass(outDat) <- "Surve"
    invisible(outDat)
}


#' An S4 Class to Represent a Fitted Model
#'
#' \code{survEm-class} is an S4 class that represents a fitted model.
#' Function \code{\link{survEm}} produces objects of this class.
#' See ``Slots'' for details.
#'
#' @slot call Function call.
#' @slot formula Formula.
#' @slot data A data frame.
#' @slot nObs A positive integer.
#' @slot estimates A list.
#' @slot control A list.
#' @slot start A list.
#' @slot na.action A length-one character vector.
#' @slot xlevels A list.
#' @slot contrasts A list.
#' @slot convergCode A nonnegative integer.
#' @slot partLogL A numeric value.
#' @slot logL A numeric value.
#' @slot fisher A numeric matrix.
#' @aliases survEm-class
#' @seealso \code{\link{survEm}} for details of slots.
#' @export
setClass(Class = "coxEm",
         slots = c(call = "call",
                   formula = "formula",
                   data = "data.frame",
                   nObs = "integer",
                   estimates = "list",
                   control = "list",
                   start = "list",
                   na.action = "character",
                   xlevels = "list",
                   contrasts = "list",
                   convergCode = "integer",
                   partLogL = "numeric",
                   logL = "numeric",
                   fisher = "matrix"))
