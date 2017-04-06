##' An S4 Class Representing Formula Response for Integrative Survival Model
##'
##' @slot id Identificator of each subject.
##' @slot time Time of reccurence event or censoring.
##' @slot event The status indicator, 0 = censored, 1 = event.
##' @aliases survi-class
##' @export
setClass("survi", contains = "matrix",
         slots = c(id = "character",
                   time = "numeric",
                   event = "integer"))


##' An S4 Class to Represent a Fitted Model
##'
##' \code{intCox-class} is an S4 class that represents a fitted model.
##' Function \code{\link{intCox}} produces objects of this class.
##' See ``Slots'' for details.
##'
##' @slot call Function call.
##' @slot formula Formula.
##' @slot data A data frame.
##' @slot nObs A positive integer.
##' @slot estimates A list.
##' @slot control A list.
##' @slot start A list.
##' @slot na.action A length-one character vector.
##' @slot xlevels A list.
##' @slot contrasts A list.
##' @slot convergCode A nonnegative integer.
##' @slot partLogL A numeric value.
##' @slot logL A numeric value.
##' @slot fisher A numeric matrix.
##' @aliases intCox-class
##' @seealso \code{\link{intCox}} for details of slots.
##' @export
setClass(Class = "intCox",
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


##' An S4 Class to Represent a Summarized Model
##'
##' \code{summaryCoxEm-class} is an S4 class that represents a summarized model.
##' Function \code{\link{summary_intCox}} produces objects of this class.  See
##' ``Slots'' for details.
##'
##' @slot call Function call.
##' @slot coefMat A matrix.
##' @aliases summaryCoxEm-class
##' @seealso \code{\link{summaryCoxEm}} for details of slots.
##' @export
setClass(Class = "summary_intCox",
         slots = c(call = "call",
                   coefMat = "matrix"))
