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


##' Concordance Index
##'
##' Compute concordance index (C-index or C-statistic) that allows weights for
##' right-censored survival data.  Asano and Hirakawa (2017) proposed using the
##' probability of being susceptible as weights for cure models.  It reduces to
##' Harrell's C-index for equal weights.
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
##' @param event A numeric vector for event indicators.  If it is \code{NULL}
##'     (by default) or \code{NA}, \code{event} will be treated all as ones and
##'     the function will compute concordance index for uncensored survival
##'     data.
##' @param risk_score A numeric vector for risk scores.
##' @param weight A optional numeric vector for weights.  If it is \code{NULL}
##'     (by default) or \code{NA}, equal weights will be used.
##'
##' @return
##' A named numeric vector that consists of
##' \itemize{
##'   \item \code{index}: the concordance index.
##'   \item \code{concordant}: the number of concordant pairs.
##'   \item \code{comparable}: the number of comparable pairs.
##'   \item \code{tied_tisk}: the number of comparable pairs having tied risks.
##' }
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
##' @examples
##' ## See examples of function 'cox_cure'.
##' @export
cIndex <- function(time, event = NULL, risk_score, weight = NULL)
{
    if (is.null(weight) || is.na(weight)) {
        weight <- 1
    }
    if (is.null(event) || is.na(event)) {
        event <- 1
    } else if (! any(event > 0)) {
        stop("No compariable pairs can be found")
    }
    rcpp_cIndex(time, event, risk_score, weight)
}


##' Bayesian Information Criterion (BIC)
##'
##' Compute Bayesian information criterion (BIC) or Schwarz's Bayesian criterion
##' (SBC) for possibly one or several objects.
##'
##' @param object An object for a fitted model.
##' @param method A character string specifying the method for computing the BIC
##'     values.  The available options are \code{"obs"} for regular BIC based on
##'     the number of observations, and \code{"effective"} for using BIC based
##'     on the number of effective sample size for censored data (number of
##'     uncensored events) proposed by Volinsky and Raftery (2000).  The former
##'     is used by default.
##' @param ... Other objects.
##'
##' @references
##'
##' Volinsky, C. T., & Raftery, A. E. (2000). Bayesian information criterion for
##' censored survival models. Biometrics, 56(1), 256--262.
##'
##' @examples
##' ## See examples of function 'cox_cure'.
##' @importFrom stats BIC
##' @export
BIC.cox_cure <- function(object, ..., method = c("obs", "effective"))
{
    method <- match.arg(method)
    bic_name <- switch(method, "obs" = "bic1", "effective" = "bic2")
    if (! missing(...)) {
        inpList <- list(object, ...)
        ## check on object class
        checkRes <- sapply(inpList, is_cox_cure)
        if (any(! checkRes))
            stop("All objects must be of the 'cox_cure' class.")
        bics <- sapply(inpList, function(a) a$model[[bic_name]])
        dfs <- sapply(inpList, function(a) a$model$coef_df)
        val <- data.frame(df = dfs, BIC = bics)
        Call <- match.call()
        is_obj <- names(Call) != "method"
        row.names(val) <- as.character(Call[is_obj][- 1L])
        return(val)
    }
    ## else return
    object$model[[bic_name]]
}


##' Bayesian Information Criterion (BIC)
##'
##' Compute Bayesian information criterion (BIC) or Schwarz's Bayesian criterion
##' (SBC) for possibly one or several objects.
##'
##' @param object An object for a fitted model.
##' @param method A character string specifying the method for computing the BIC
##'     values.  The available options are \code{"obs"} for regular BIC based on
##'     the number of observations, and \code{"certain-event"} for a variant of
##'     BIC based on the number of certain uncensored events.  The former is
##'     used by default.
##' @param ... Other objects.
##'
##' @examples
##' ## See examples of function 'cox_cure'.
##' @importFrom stats BIC
##' @export
BIC.cox_cure_uncer <- function(object, ..., method = c("obs", "certain-event"))
{
    method <- match.arg(method)
    bic_name <- switch(method, "obs" = "bic1", "certain-event" = "bic2")
    if (! missing(...)) {
        inpList <- list(object, ...)
        ## check on object class
        checkRes <- sapply(inpList, is_cox_cure_uncer)
        if (any(! checkRes))
            stop("All objects must be of the 'cox_cure_uncer' class.")
        bics <- sapply(inpList, function(a) a$model[[bic_name]])
        dfs <- sapply(inpList, function(a) a$model$coef_df)
        val <- data.frame(df = dfs, BIC = bics)
        Call <- match.call()
        is_obj <- names(Call) != "method"
        row.names(val) <- as.character(Call[is_obj][- 1L])
        return(val)
    }
    ## else return
    object$model[[bic_name]]
}


##' Bayesian Information Criterion (BIC)
##'
##' Compute Bayesian information criterion (BIC) or Schwarz's Bayesian criterion
##' (SBC) from a fitted model.
##'
##' @param object An object for a fitted model.
##' @param ... Other arguments for future usage.
##'
##' @examples
##' ## See examples of function 'cox_cure_net'.
##' @importFrom stats BIC
##' @export
BIC.cox_cure_net <- function(object, ...)
{
    bics <- object$model$bic
    dfs <- object$model$coef_df
    data.frame(df = dfs, BIC = bics)
}
