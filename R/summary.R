##
## intsurv: Integrative Survival Models
## Copyright (C) 2017-2021  Wenjie Wang <wang@wwenjie.org>
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


##' Summary of a Fitted Model
##'
##' For \code{\link{iCoxph}} object, the function returns a
##' \code{\link{summary.iCoxph}} object whose slots include
##' \itemize{
##'     \item \code{call}: Function call of model fitting.
##'     \item \code{coefMat}: Estimated covariate coefficients.
##'     \item \code{logL}: Log-likelihood under observed data.
##' }
##'
##' @param object \code{\link{iCoxph-class}} object.
##' @param showCall A logic value with default \code{TRUE}, indicating whether
##'     \code{show} method prints out the original call information of
##'     \code{iCoxph}.  Set \code{FALSE} for a more concise printout.
##' @param ... Other arguments for future usage.
##'
##' @return \code{summary.iCoxph} class object.
##'
##' @aliases summary.iCoxph summary,iCoxph-method
##'
##' @examples
##' ## See examples of function iCoxph
##' @seealso
##' \code{\link{iCoxph}} for model fitting;
##' \code{\link{coef,iCoxph-method}} for coefficient estimates.
##' @export
setMethod(
    f = "summary", signature = "iCoxph",
    definition = function(object, showCall = TRUE, ...)
    {
        Call <- object@call
        attr(Call, "show") <- showCall
        betaMat <- object@estimates$beta
        len_logL <- length(object@logL)
        new("summary.iCoxph",
            call = Call,
            coefMat = betaMat,
            logL = object@logL[len_logL])
    }
)


##' Summary of a Fitted Model
##'
##' Summarize a fitted Cox cure rate model with possible uncertain event status.
##'
##' @param object Object representing a fitted model.
##' @param ... Other arguments for future usage.  A warning will be thrown if
##'     any invalid argument is specified.
##'
##' @method summary cox_cure
##' @importFrom stats pnorm
##' @export
summary.cox_cure <- function(object, ...)
{
    warn_dots()
    ## function summarizing coef matrix
    get_coef_mat <- function(object, part = c("surv", "cure"), ...) {
        part <- match.arg(part)
        coef_name <- sprintf("%s_coef", part)
        se_name <- sprintf("%s_coef_se", part)
        mat <- matrix(NA, nrow = length(object[[coef_name]]), ncol = 5L)
        colnames(mat) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
        rownames(mat) <- names(object[[coef_name]])
        mat[, "coef"] <- object[[coef_name]]
        mat[, "exp(coef)"] <- exp(mat[, "coef"])
        mat[, "se(coef)"] <- if (is.null(object$bootstrap[[se_name]])) {
                                 NA_real_
                             } else {
                                 object$bootstrap[[se_name]]
                             }
        mat[, "z"] <- mat[, "coef"] / mat[, "se(coef)"]
        mat[, "Pr(>|z|)"] <- 2 * stats::pnorm(- abs(mat[, "z"]))
        ## return
        mat
    }
    ## return
    structure(
        list(
            surv_coef_mat = get_coef_mat(object, "surv"),
            cure_coef_mat = get_coef_mat(object, "cure"),
            model = object$model,
            call = object$call
        ),
        class = "summary_cox_cure"
    )
}


##' @rdname summary.cox_cure
##' @method summary cox_cure_mar
##' @export
summary.cox_cure_mar <- function(object, ...)
{
    warn_dots()
    ## function summarizing coef matrix
    get_coef_mat <- function(object, part = c("surv", "cure"), ...) {
        part <- match.arg(part)
        coef_name <- sprintf("%s_coef", part)
        se_name <- sprintf("%s_coef_se", part)
        mat <- matrix(NA, nrow = length(object[[coef_name]]), ncol = 5L)
        colnames(mat) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
        rownames(mat) <- names(object[[coef_name]])
        mat[, "coef"] <- object[[coef_name]]
        mat[, "exp(coef)"] <- exp(mat[, "coef"])
        mat[, "se(coef)"] <- if (is.null(object$bootstrap[[se_name]])) {
                                 NA_real_
                             } else {
                                 object$bootstrap[[se_name]]
                             }
        mat[, "z"] <- mat[, "coef"] / mat[, "se(coef)"]
        mat[, "Pr(>|z|)"] <- 2 * stats::pnorm(- abs(mat[, "z"]))
        ## return
        mat
    }
    ## return
    structure(
        list(
            surv_coef_mat = get_coef_mat(object, "surv"),
            cure_coef_mat = get_coef_mat(object, "cure"),
            model = object$model,
            call = object$call
        ),
        class = "summary_cox_cure_mar"
    )
}
