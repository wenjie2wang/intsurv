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


##' @method print cox_cure
##' @export
print.cox_cure <- function(x, ...) {
    ## function call
    cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")
    ## survival part
    surv_mat <- matrix(NA, nrow = length(x$surv_coef), ncol = 5L)
    colnames(surv_mat) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
    rownames(surv_mat) <- names(x$surv_coef)
    surv_mat[, "coef"] <- x$surv_coef
    surv_mat[, "exp(coef)"] <- exp(surv_mat[, "coef"])
    surv_mat[, "se(coef)"] <-
        if (is.null(x$bootstrap$surv_coef_se)) {
            NA_real_
        } else {
            x$bootstrap$surv_coef_se
        }
    surv_mat[, "z"] <- surv_mat[, "coef"] / surv_mat[, "se(coef)"]
    surv_mat[, "Pr(>|z|)"] <- 2 * stats::pnorm(- abs(surv_mat[, "z"]))
    cat("\nCoefficient estimates from the survival part:\n\n")
    printCoefmat(surv_mat)
    ## cure rate part
    cure_mat <- matrix(NA, nrow = length(x$cure_coef), ncol = 5L)
    colnames(cure_mat) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
    rownames(cure_mat) <- names(x$cure_coef)
    cure_mat[, "coef"] <- x$cure_coef
    cure_mat[, "exp(coef)"] <- exp(cure_mat[, "coef"])
    cure_mat[, "se(coef)"] <-
        if (is.null(x$bootstrap$cure_coef_se)) {
            NA_real_
        } else {
            x$bootstrap$cure_coef_se
        }
    cure_mat[, "z"] <- cure_mat[, "coef"] / cure_mat[, "se(coef)"]
    cure_mat[, "Pr(>|z|)"] <- 2 * stats::pnorm(- abs(cure_mat[, "z"]))
    cat("\nCoefficient estimates from the cure rate part:\n\n")
    printCoefmat(cure_mat)
    ## return
    invisible(x)
}

## for cox_cure_uncer
print.cox_cure_uncer <- print.cox_cure
