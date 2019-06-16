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

## collation after class.R
##' @include class.R
NULL


##' Weighted Concordance Index
##'
##' This function is a straightforward implementation of concordance index
##' (C-index or C-statistic) that allows weights.  Asano and Hirakawa (2017)
##' proposed using the probability of being susceptible as weights for cure
##' models.  It reduces to Harrell's C-index for equal weights.
##'
##' Let \eqn{r_i}, \eqn{t_i}, and \eqn{\delta_i} denote the risk score, observed
##' time, and event indicator of \eqn{i}-th subject.  The pair of
##' \eqn{(t_i,\delta_i)} and \eqn{(t_j,\delta_j)}, where \eqn{i<j}, are defined
##' to be comparable if \eqn{t_i<t_j,\delta_i=1} or
##' \eqn{t_i=t_j,\delta_i=1,\delta_j=0}.  In the context of survival analysis,
##' the risk scores of a comparable pair are said to be concordant with the
##' survival outcomes if \eqn{r_i>r_j}.  The C-index is defined as the
##' proportion of the concordant pairs among the comparable pairs.  For
##' comparable pair satisfying \eqn{t_i<t_j,\delta_i=1}, we count 0.5 in the
##' numerator of the concordance index for tied risk scores (\eqn{r_i=r_j}).
##'
##' @param time A numeric vector for observed times
##' @param event A numeric vector for event indicators.
##' @param risk_score A numeric vector for risk scores.
##' @param weight A optional numeric vector for weights.  By default, equal
##'     weights are used.
##'
##' @references
##'
##' Asano, J., & Hirakawa, A. (2017). Assessing the prediction accuracy of a
##' cure model for censored survival data with long-term survivors: Application
##' to breast cancer data. Journal of biopharmaceutical statistics, 27(6),
##' 918--932.
##'
##' Harrell, F. E., Lee, K. L., & Mark, D. B. (1996). Multivariable prognostic
##' models: Issues in developing models, evaluating assumptions and adequacy,
##' and measuring and reducing errors. Statistics in medicine, 15(4), 361--387.
##'
##' @export
cIndex <- function(time, event, risk_score, weight = NULL) {
    if (is.null(weight) || is.na(weight)) {
        weight <- 1
    }
    rcpp_cIndex(time, event, risk_score, weight)
}
