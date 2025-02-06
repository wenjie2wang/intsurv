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


##' @method print cox_cure
##' @export
print.cox_cure <- function(x, ...)
{
    if (! is.null(x$call)) {
        ## function call
        cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n", sep = "")
    }
    ## get the coef matrix
    obj <- summary.cox_cure(x)
    ## survival part
    cat("\nCoefficient estimates for the survival part:\n\n")
    printCoefmat(obj$surv_coef_mat)
    ## cure rate part
    cat("\nCoefficient estimates for the cure rate part:\n\n")
    printCoefmat(obj$cure_coef_mat)
    ## return
    invisible(x)
}


##' @method print cox_cure_mar
##' @export
print.cox_cure_mar <- function(x, ...)
{
    if (! is.null(x$call)) {
        ## function call
        cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n", sep = "")
    }
    ## get the coef matrix
    obj <- summary.cox_cure_mar(x)
    ## survival part
    cat("\nCoefficient estimates for the survival part:\n\n")
    printCoefmat(obj$surv_coef_mat)
    ## cure rate part
    cat("\nCoefficient estimates for the cure rate part:\n\n")
    printCoefmat(obj$cure_coef_mat)
    ## return
    invisible(x)
}


##' @method print summary_cox_cure
##' @export
print.summary_cox_cure <- function(x, ...)
{
    if (! is.null(x$call)) {
        ## function call
        cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n", sep = "")
    }
    ## survival coef
    cat("\nCoefficient estimates for the survival part:\n\n")
    printCoefmat(x$surv_coef_mat)
    ## cure coef
    cat("\nCoefficient estimates for the cure rate part:\n\n")
    printCoefmat(x$cure_coef_mat)
    ## misc
    cat(sprintf("\nNumber of observations: %d\n", x$model$nObs))
    cat(sprintf("Number of events: %d\n", x$model$nEvent))
    cat(sprintf("Weighted concordance index: %2.2f%s\n",
                100 * x$model$c_index, "%"))
    cat(sprintf("Observed data log-likelihood: %.3f\n",
                - x$model$negLogL))
    ## return
    invisible(x)
}


## TODO: customize more for cox_cure_mar
##' @method print summary_cox_cure_mar
##' @export
print.summary_cox_cure_mar <- print.summary_cox_cure


##' @method print cox_cure_net
##' @export
print.cox_cure_net <- function(x, ...)
{
    if (! is.null(x$call)) {
        ## function call
        cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
    }
    bic_dat <- BIC.cox_cure_net(x)
    digits <- max(3, getOption("digits") - 3)
    lambda_dat <- as.data.frame(signif(x$penalty$lambda_mat, digits))
    colnames(lambda_dat) <- c("surv_l1", "surv_l2", "cure_l1", "cure_l2")
    print_dat <- cbind(bic_dat, lambda_dat)
    print(print_dat)
    ## return nothing
    invisible(NULL)
}


##' @method print cox_cure_net_mar
##' @export
print.cox_cure_net_mar <- print.cox_cure_net
