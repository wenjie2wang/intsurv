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
##' bootSe(object, B = 50, se = c("mad", "inter-quantile", "sd"), ...)
##'
##' @param object \code{\link{iCoxph-class}} object.
##' @param B A positive integer specifying number of bootstrap samples
##'     used for SE estimates.  A large number, such as 200, is often needed for
##'     a more reliable estimation in practice.
##' @param se A character value specifying the way computing SE from bootstrap
##'     samples. The default method is based on median absolute deviation and
##'     the second method is based on inter-quantile, both of which are based on
##'     normality of the bootstrap estimates and provids robust estimates for
##'     SE. The third method estimates SE by the standard deviation of the
##'     bootstrap estimates.
##' @param ... Other arguments for future usage.
##'
##' @return \code{\link{iCoxph-class}} object.
##'
##' @examples
##' ## See examples given in function 'iCoxph'
##' @seealso
##' \code{\link{iCoxph}} for fitting integerative Cox model.
##' @importFrom stats median pnorm qnorm quantile sd
##' @export
bootSe <- function(object,
                   B = 50,
                   se = c("mad", "inter-quantile", "sd"),
                   ...)
{
    if (! is_iCoxph(object))
        stop("The 'object' has to be an 'iCoxph' class object.")
    if (! is.integer(B))
        B <- as.integer(B)
    if (B <= 0)
        stop("The number of bootstrap samples must be a postive integer.")
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
    ## if (control$estOnly) {
    ##     nBeta <- NROW(estMat)
    ##     estMat <- t(estMat)
    ##     colnames(estMat) <- paste0("b", seq_len(nBeta))
    ##     return(estMat)
    ## }
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
    object@estimates$boostrap_se <- cbind(
        "mad" = se_mad,
        "inter-quantile" = se_interQ,
        "sd" = se_sd
    )
    tmp <- object@estimates$beta[, "z"] <-
        object@estimates$beta[, "coef"] / object@estimates$beta[, "se(coef)"]
    object@estimates$beta[, "Pr(>|z|)"] <- 2 * stats::pnorm(- abs(tmp))

    ## return
    object
}
