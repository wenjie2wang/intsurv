##
## intsurv: Integrative Survival Models
## Copyright (C) 2017-2018  Wenjie Wang <wjwang.stat@gmail.com>
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
    cal$start <- start
    ## update local control list
    control <- do.call(bootSe_control, control)
    ## add noSE = TRUE to the original control list
    fm <- object@formula
    cal$control <- object@control
    cal$control$noSE <- TRUE
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
    estMat <- replicate(numBoot, {
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
    if (control$estOnly) {
        nBeta <- NROW(estMat)
        estMat <- t(estMat)
        colnames(estMat) <- paste0("b", seq_len(nBeta))
        return(estMat)
    }
    ## some computing can be skipped but kept now for testing
    se_mad <- apply(estMat, 1L, function(a) {
            median(abs(a - median(a))) * 1.4826
    })
    se_interQ <- apply(estMat, 1L, function(a) {
            diff(stats::quantile(a, probs = c(0.25, 0.75))) /
                (stats::qnorm(0.75) - stats::qnorm(0.25))
    })
    se_sd <- apply(estMat, 1L, sd)
    object@estimates$beta[, "se(coef)"] <-
        switch(se,
               "mad" = se_mad,
               "inter-quantile" = se_interQ,
               "sd" = se_sd)

    ## save estMat for testing
    object@estimates$boostrap_beta <- estMat
    object@estimates$boostrap_se <- cbind("mad" = se_mad,
                                          "inter-quantile" = se_interQ,
                                          "sd" = se_sd)

    tmp <- object@estimates$beta[, "z"] <- object@estimates$beta[, "coef"] /
        object@estimates$beta[, "se(coef)"]
    object@estimates$beta[, "Pr(>|z|)"] <- 2 * stats::pnorm(- abs(tmp))
    object
}


### internal functions =========================================================
bootSe_start <- function(betaVec = NULL,
                         betaMat = NULL,
                         piVec = NULL,
                         censorRate = NULL,
                         ...,
                         start0)
{
    if (is.null(betaVec))
        betaVec <- start0$beta0
    if (is.null(piVec))
        piVec <- start0$piVec

    censorRate0 <- start0$censorRate0
    if (is.null(censorRate) && ! is.na(censorRate0)) {
        censorRate <- if (start0$multiStart) {
                          seq.int(max(0, censorRate0 - 0.2),
                                  min(1, censorRate0 + 0.2), 0.05)
                      } else {
                          start0$censorRate
                      }
    }

    ## return
    list(betaVec = betaVec,
         betaMat = betaMat,
         piVec = piVec,
         censorRate = censorRate,
         semiparametric = start0$semiparametric,
         parametric = start0$parametric,
         parametricOnly = start0$parametricOnly,
         multiStart = start0$multiStart,
         randomly = start0$randomly)
}


bootSe_control <- function(estOnly = FALSE, ...)
{
    if (! is.logical(estOnly))
        stop("Argument 'estOnly' has to be a logical value.")
    list(estOnly = estOnly)
}
