################################################################################
##
##   R package intsurv by Wenjie Wang, Kun Chen and Jun Yan
##   Copyright (C) 2017-2018
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


## collation after class.R
##' @include class.R
NULL


##' Standard Error Estimates through Bootstrapping Methods
##'
##' This function addes or updates standard error (SE) estimates through
##' bootstrap method for \code{\link{iCoxph-class}} object by default.  Three
##' different methods are available for computing SE from bootstrap samples
##' through argument \code{se}.
##'
##' Given the fact that the bootstrap method is computional intensive, the
##' function can return the coefficient estimates from one bootstrap sample when
##' \code{control = list(estOnly = TRUE)} is specified, which can used in
##' parallel computing or high performance computing (HPC) cluster. Then the SE
##' estimates can be easily computed based on estimates from bootstrap samples.
##'
##' The available elements of argument \code{start} are the same with those of
##' argument \code{start} in function \code{\link{iCoxph}} except that the
##' \code{piVec} is not available since its length may vary in different
##' bootstrap samples.
##'
##' @usage
##' bootSe(object, numBoot = 50, se = c("mad", "inter-quantile", "sd"),
##'        start = list(), control = list(), ...)
##'
##' @param object \code{\link{iCoxph-class}} object.
##' @param numBoot A positive integer specifying number of bootstrap samples
##'     used for SE estimates.  A large number, such as 200, is often needed for
##'     a more reliable estimation in practice.
##' @param se A character value specifying the way computing SE from bootstrap
##'     samples. The default method is based on median absolute deviation and
##'     the second method is based on inter-quantile, both of which are based on
##'     normality of the bootstrap estimates and provids robust estimates for
##'     SE. The third method estimates SE by the standard deviation of the
##'     bootstrap estimates.
##' @param start An optional list of starting values for the parameters to be
##'     estimated in the model.  See more in Section details.
##' @param control A optional list that controls the model fitting process for
##'     bootstrap samples, or specifys the function to return coefficient
##'     estimates from the bootstrap samples. See the available options in
##'     Section Details.
##' @param ... Other arguments for future usage.
##'
##' @return \code{\link{iCoxph-class}} object by default or a numeric matrix of
##'     coefficient estimates from each bootstrap sample.
##'
##' @examples
##' ## See examples given in function 'iCoxph'
##' @seealso
##' \code{\link{iCoxph}} for fitting extended Cox model for uncertain records.
##' @importFrom stats median pnorm qnorm quantile sd
##' @export
bootSe <- function(object, numBoot = 50, se = c("mad", "inter-quantile", "sd"),
                   start = list(), control = list(), ...)
{
    if (! inherits(object, "iCoxph"))
        stop("The 'object' has to be an 'iCoxph' class object.")
    se <- match.arg(se)
    cal <- object@call
    ## update start list
    start <- do.call(bootSe_start, c(start, list(start0 = object@start)))
    cal$start <- quote(start)
    ## update local control list
    control <- do.call(bootSe_control, control)
    ## add noSE = TRUE to the original control list
    fm <- object@formula
    cal$control <- object@control
    cal$control$noSE <- TRUE
    cal$data <- quote(bootDat)
    dat <- object@data
    idName <- as.character(fm[[2L]][[2L]])
    timeName <- as.character(fm[[2L]][[3L]])
    eventName <- as.character(fm[[2L]][[4L]])
    dat <- dat[order(dat[, idName], dat[, timeName],
                     1L - dat[, eventName]), ]
    uid <- unique(as.character(dat[, idName]))
    idTab <- table(dat[, idName])
    estMat <- replicate(numBoot, {
        sID <- sort(sample(uid, replace = TRUE))
        tmpDat <- data.frame(ID = sID)
        colnames(tmpDat) <- idName
        repNum <- idTab[match(sID, uid)]
        bootDat <- merge(tmpDat, dat, by = idName)
        bootDat[[idName]] <- rep(seq_along(sID), times = repNum)
        res <- eval(cal)
        as.numeric(res@estimates$beta[, "coef"])
    })
    if (control$estOnly) {
        nBeta <- NROW(estMat)
        estMat <- t(estMat)
        colnames(estMat) <- paste0("b", seq_len(nBeta))
        return(estMat)
    }
    if (identical(se, "mad")) {
        object@estimates$beta[, "se(coef)"] <- apply(estMat, 1L, function(a) {
            median(abs(a - median(a))) * 1.4826
        })
    } else if (identical(se, "inter-quantile")) {
        object@estimates$beta[, "se(coef)"] <- apply(estMat, 1L, function(a) {
            diff(stats::quantile(a, probs = c(0.25, 0.75))) /
                (stats::qnorm(0.75) - stats::qnorm(0.25))
        })
    } else (identical(se, "sd"))
        object@estimates$beta[, "se(coef)"] <- apply(estMat, 1L, sd)

    tmp <- object@estimates$beta[, "z"] <- object@estimates$beta[, "coef"] /
        object@estimates$beta[, "se(coef)"]
    object@estimates$beta[, "Pr(>|z|)"] <- 2 * stats::pnorm(- abs(tmp))
    object
}


### internal functions =========================================================
bootSe_start <- function(beta0, censorRate, multiStart = TRUE, ...,
                         start0)
{
    censorRate0 <- start0$censorRate0
    if (missing(beta0))
        beta <- start0$beta
    if (missing(censorRate)) {
        censorRate <- if (! multiStart) {
                          censorRate0
                      } else {
                          seq.int(max(0, censorRate0 - 0.2),
                                  min(1, censorRate0 + 0.2), 0.05)
                      }
    }
    list(beta = beta, censorRate = censorRate,
         censorRate0 = censorRate0)
}


bootSe_control <- function(estOnly = FALSE, ...)
{
    if (! is.logical(estOnly))
        stop("Argument 'estOnly' has to be a logical value.")
    list(estOnly = estOnly)
}
