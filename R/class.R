##
## intsurv: Integrative Survival Models
## Copyright (C) 2017-2025  Wenjie Wang <wang@wwenjie.org>
##
## This file is part of the R package intsurv.
##
## The R package intsurv is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package intsurv is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

##' An S4 Class Representing Formula Response
##'
##' @slot .Data A numeric matrix object.
##' @slot ID Identificator of each subject.
##' @aliases Survi-class
##' @export
setClass("Survi", contains = "matrix",
         slots = c(ID = "character"))


##' An S4 Class to Represent a Fitted Integrative Cox Model
##'
##' The \code{iCoxph} class is an S4 class that represents a fitted model.
##' Function \code{\link{iCoxph}} produces objects of this class.  See ``Slots''
##' for details.
##'
##' @slot call Function call.
##' @slot formula Model formula.
##' @slot nObs A positive integer.
##' @slot data A data frame.
##' @slot estimates A list.
##' @slot start A list.
##' @slot control A list.
##' @slot na.action A length-one character vector.
##' @slot xlevels A list.
##' @slot contrasts A list.
##' @slot convergeCode A non-negative integer.
##' @slot logL A numeric value.
##'
##' @aliases iCoxph-class
##'
##' @seealso \code{\link{iCoxph}} for details of slots.
##' @export
setClass(Class = "iCoxph",
         slots = c(call = "call",
                   formula = "formula",
                   data = "data.frame",
                   nObs = "integer",
                   estimates = "list",
                   start = "list",
                   control = "list",
                   na.action = "character",
                   xlevels = "list",
                   contrasts = "list",
                   convergeCode = "integer",
                   logL = "numeric"))


##' An S4 Class to Represent Summary of a fitted integrative Cox model
##'
##' The \code{summary.intsurv} class is an S4 class that represents a summarized
##' model.  Function \code{\link{summary,iCoxph-method}} produces objects of
##' this class.  See ``Slots'' for details.
##'
##' @slot call Function call.
##' @slot coefMat A matrix.
##' @slot logL A numeric value.
##'
##' @aliases summary.iCoxph-class
##'
##' @seealso \code{\link{summary,iCoxph-method}} for meaning of slots.
##'
##' @export
setClass(Class = "summary.iCoxph",
         slots = c(call = "call",
                   coefMat = "matrix",
                   logL = "numeric"))
