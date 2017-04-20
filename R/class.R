################################################################################
##
##   R package intsurv by Wenjie Wang, Kun Chen and Jun Yan
##   Copyright (C) 2017
##
##   This file is part of the R package intsurv.
##
##   The R package intsurv is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package intsurv is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################


##' An S4 Class Representing Formula Response
##'
##' @slot .Data A numeric matrix object.
##' @slot id Identificator of each subject.
##' @slot time Time of reccurence event or censoring.
##' @slot event The status indicator, 0 = censored, 1 = event.
##' @aliases Survi-class
##' @export
setClass("Survi", contains = "matrix",
         slots = c(id = "character",
                   time = "numeric",
                   event = "integer"))


##' An S4 Class to Represent a Fitted Extended Cox Model
##'
##' \code{coxphx-class} is an S4 class that represents a fitted model.
##' Function \code{\link{coxphx}} produces objects of this class.
##' See ``Slots'' for details.
##'
##' @slot call Function call.
##' @slot formula Formula.
##' @slot nObs A positive integer.
##' @slot data A data frame.
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
##' @aliases coxphx-class
##' @seealso \code{\link{coxphx}} for details of slots.
##' @export
setClass(Class = "coxphx",
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


##' An S4 Class to Represent Summary of a Fitted Model
##'
##' \code{coxphx.summary-class} is an S4 class that represents a summarized
##' model.  Function \code{\link{summary,coxphx-method}} produces objects of
##' this class.  See ``Slots'' for details.
##'
##' @slot call Function call.
##' @slot coefMat A matrix.
##' @slot logL A numeric value.
##' @aliases coxphx.summary-class
##' @seealso \code{\link{summary,coxphx-method}} for meaning of slots.
##' @export
setClass(Class = "coxphx.summary",
         slots = c(call = "call",
                   coefMat = "matrix",
                   logL = "numeric"))
