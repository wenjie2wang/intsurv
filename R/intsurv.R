##
## intsurv: Integrative Survival Models
## Copyright (C) 2017-2020  Wenjie Wang <wjwang.stat@gmail.com>
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
##' \item integrative Cox model with uncertain event times (Wang et al., 2020)
##'
##' \item Cox cure rate model with uncertain event status (Wang et al., 2020)
##'
##' }
##'
##' It also contains other survival analysis routines, including regular Cox
##' cure rate model, regularized Cox cure rate model with elastic net penalty,
##' and weighted concordance index.
##'
##' @references
##'
##' Wang, W., Aseltine, R. H., Chen, K., & Yan, J. (2020). Integrative Survival
##' Analysis with Uncertain Event Times in Application to A Suicide Risk
##' Study. \emph{Annals of Applied Statistics}, 14(1), 51--73.
##'
##' Wang, W., Luo, C., Aseltine, R. H., Wang, F., Yan, J., & Chen, K. (2020).
##' Suicide Risk Modeling with Uncertain Diagnostic Records. \emph{arXiv
##' preprint arXiv:2009.02597}.
##'
##' @importFrom Rcpp sourceCpp
##' @useDynLib intsurv
##'
##' @importFrom methods setClass setGeneric setMethod new validObject
##'
##' @docType package
##' @name intsurv-package
NULL
