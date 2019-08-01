##
## intsurv: Integrative Survival Models
## Copyright (C) 2017-2019  Wenjie Wang <wjwang.stat@gmail.com>
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


##' Integrative Survival Modeling
##'
##' The package \pkg{intsurv} provides implementations of
##'
##' \itemize{
##'
##' \item integrative Cox model with uncertain event times (Wang et al., 2019)
##'
##' \item Cox cure rate model with uncertain event status (Wang et al., 2019+)
##'
##' }
##'
##' It also contains other survival analysis routines, including regular Cox
##' cure rate model, regularized Cox cure rate model with elastic net penalty,
##' and weighted concordance index.
##'
##' @references
##'
##' Wang, W., Aseltine, R., Chen, K., & Yan, J. (2019).  Integrative Survival
##' Analysis with Uncertain Event Times in Application to a Suicide Risk
##' Study. \emph{Annals of Applied Statistics}. (in press)
##'
##' Wang, W., Chen, K., Luo, C., & Yan, J. (2019+). Cox Cure Model with
##' Uncertain Event Status with application to a Suicide Risk
##' Study. \emph{Working in Progress}.
##'
##' @importFrom Rcpp sourceCpp
##' @useDynLib intsurv
##'
##' @importFrom methods setClass setGeneric setMethod new validObject
##'
##' @docType package
##' @name intsurv-package
NULL
