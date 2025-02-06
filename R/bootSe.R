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

## collation after class.R
##' @include class.R
NULL


##' Standard Error Estimates through Bootstrap
##'
##' For \code{\link{iCoxph-class}} object, add (or update) standard error (SE)
##' estimates through bootstrap methods, or compute the coefficient estimates
##' from the given number of bootstrap samples.
##'
##' Three different methods are available for computing SE from bootstrap
##' samples through argument \code{se}.  Given the fact that the bootstrap
##' method is computationally intensive, the function returns the coefficient
##' estimates in a matrix from the given number of bootstrap samples when
##' \code{return_beta = TRUE)} is specified, which can be used in parallel
##' computing or high performance computing (HPC) cluster.  The SE estimates can
##' be further computed based on estimates from bootstrap samples by users on
##' their own.  The \code{return_beta = TRUE} is implied, when \code{B = 1} is
##' specified.
##'
##' @param object \code{\link{iCoxph-class}} object.
##' @param B A positive integer specifying number of bootstrap samples used for
##'     SE estimates.  A large number, such as 200, is often needed for a more
##'     reliable estimation in practice.  If \code{B = 1} is specified, the
##'     function will return the covariate coefficient estimates instead of a
##'     \code{iCoxph-class} object.
##' @param se A character value specifying the way computing SE from bootstrap
##'     samples. The default method is based on median absolute deviation and
##'     the second method is based on inter-quartile, both of which are based on
##'     normality of the bootstrap estimates and provides robust estimates for
##'     SE. The third method estimates SE by the standard deviation of the
##'     bootstrap estimates.
##' @param return_beta A logical value. If \code{TRUE}, the function returns the
##'     covariate coefficient estimates from the given number of bootstrap
##'     samples, which allows users to split the most computationally intensive
##'     step into small pieces that can be computed in a parallel manner.  The
##'     default value is \code{FALSE}.
##' @param ... Other arguments for future usage.  A warning will be thrown if
##'     any invalid argument is specified.
##'
##' @return \code{\link{iCoxph-class}} object or a numeric matrix that contains
##'     the covariate coefficient estimates from the given number of bootstrap
##'     samples in rows.
##'
##' @examples
##' ## See examples of function 'iCoxph'.
##' @seealso
##' \code{\link{iCoxph}} for fitting integrative Cox model.
##' @importFrom stats median pnorm qnorm quantile sd
##' @export
bootSe <- function(object,
                   B = 50,
                   se = c("inter-quartile", "mad", "sd"),
                   return_beta = FALSE,
                   ...)
{
    ## some simple checkings
    if (! is_iCoxph(object))
        stop("The 'object' has to be an 'iCoxph' class object.")
    if (! is.integer(B))
        B <- as.integer(B)
    if (B <= 0)
        stop("The number of bootstrap samples must be a postive integer.")
    ## warning on `...`
    warn_dots(...)

    se <- match.arg(se)
    cal <- object@call
    fm <- object@formula
    ## pre-processing data
    cal$data <- quote(bootDat)
    dat <- object@data
    idName <- as.character(fm[[2L]][[2L]])
    timeName <- as.character(fm[[2L]][[3L]])
    eventName <- as.character(fm[[2L]][[4L]])
    dat <- dat[order(dat[, idName], dat[, timeName],
                     1L - dat[, eventName]), ]
    id_string <- as.character(dat[, idName])
    uid <- unique(id_string)
    dup_id_string <- unique(id_string[duplicated(id_string)])
    dupIdx <- id_string %in% dup_id_string
    uni_id_string <- id_string[! dupIdx]
    idTab <- table(id_string)
    estMat <- replicate(B, {
        uni_sID <- sample(uni_id_string, replace = TRUE)
        dup_sID <- sample(dup_id_string, replace = TRUE)
        sID <- sort(c(uni_sID, dup_sID))
        tmpDat <- data.frame(ID = sID, stringsAsFactors = FALSE)
        colnames(tmpDat) <- idName
        repNum <- idTab[match(sID, uid)]
        bootDat <- merge(tmpDat, dat, by = idName)
        bootDat[[idName]] <- rep(seq_along(sID), times = repNum)
        res <- eval(cal)
        as.numeric(res@estimates$beta[, "coef"])
    })

    ## return estimated beta only if B is 1
    if (B == 1L || return_beta) {
        nBeta <- NROW(estMat)
        estMat <- t(estMat)
        colnames(estMat) <- paste0("b", seq_len(nBeta))
        return(estMat)
    }

    ## compute se estimates
    se_vec <- switch(
        se,
        "inter-quartile" = {
            apply(estMat, 1L, se_interQ)
        },
        "mad" = {
            apply(estMat, 1L, function(a) {
                median(abs(a - median(a))) * 1.4826
            })
        },
        "sd" = {
            apply(estMat, 1L, sd)
        }
    )

    ## add/update se estimates in the obj
    object@estimates$beta[, "se(coef)"] <- se_vec

    ## save estMat for testing
    ## object@estimates$boostrap_beta <- estMat
    ## object@estimates$boostrap_se <- cbind(
    ##     "mad" = se_mad,
    ##     "inter-quartile" = se_interQ,
    ##     "sd" = se_sd
    ## )

    ## compute p-value based on wald test
    tmp <- object@estimates$beta[, "z"] <-
        object@estimates$beta[, "coef"] / object@estimates$beta[, "se(coef)"]
    object@estimates$beta[, "Pr(>|z|)"] <- 2 * stats::pnorm(- abs(tmp))

    ## return
    object
}
