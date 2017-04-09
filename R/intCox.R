##' Integrative Cox Model for Uncertain Survival Data
##'
##' The function fits the extended Cox model for uncertain survival data
##' due to imperfect data integration proposed by Wang (2017).
##'
##' @param formula \code{survi} object specifying the covariates and response
##'     variable in the model, such as \code{survi(ID, time, event) ~ x1 + x2}.
##' @param data An optional data frame, list, or environment that contains the
##'     covariates and response variables included in the model. If not found in
##'     data, the variables are taken from \code{environment(formula)}, usually
##'     the environment from which function \code{\link{intCox}} is called.
##' @param subset An optional vector specifying a subset of observations to be
##'     used in the fitting process.
##' @param na.action An optional function that indicates what should the
##'     procedure do if the data contains \code{NA}s.  The default is set by the
##'     na.action setting of \code{\link[base]{options}}.  The "factory-fresh"
##'     default is \code{\link[stats]{na.omit}}.  Other possible values inlcude
##'     \code{\link[stats]{na.fail}}, \code{\link[stats]{na.exclude}}, and
##'     \code{\link[stats]{na.pass}}.  \code{help(na.fail)} for details.
##' @param contrasts An optional list, whose entries are values (numeric
##'     matrices or character strings naming functions) to be used as
##'     replacement values for the contrasts replacement function and whose
##'     names are the names of columns of data containing factors.  See
##'     \code{contrasts.arg} of \code{\link[stats]{model.matrix.default}} for
##'     details.
##' @param start An optional list of starting values for the parameters to be
##'     estimated in the model.  See more in Section details.
##' @param control An optional list of parameters to control the maximization
##'     process of negative log likelihood function and adjust the baseline rate
##'     function.  See more in Section details.
##' @param ... Other arguments for future usage.
##' @return A \code{\link{intCox-class}} object, whose slots include
##' \itemize{
##'     \item \code{call}: Function call.
##'     \item \code{formula}: Formula used in the model fitting.
##'     \item \code{data}: original input dataset.
##'         (FIXME: remove this slot before submission to CRAN)
##'     \item \code{nObs}: Number of observation.
##'     \item \code{estimates}:
##'         \itemize{
##'             \item \code{spline}: The name of splines used.
##'             \item \code{knots}: Internal knots specified for the baseline
##'                 rate function.
##'             \item \code{Boundary.knots}: Boundary knots specified for the
##'                 baseline rate function.
##'             \item \code{degree}: Degree of spline bases specified in
##'                 baseline rate function.
##'             \item \code{df}: Degree of freedom of the model specified.
##'     }
##'     \item \code{estimates}: Estimated coefficients of covariates and
##'         baseline rate function, and estimated rate parameter of
##'         gamma frailty variable.
##'     \item \code{control}: The control list specified for model fitting.
##'     \item \code{start}: The initial guess specified for the parameters
##'         to be estimated.
##'     \item \code{na.action}: The procedure specified to deal with
##'         missing values in the covariate.
##'     \item \code{xlevels}: A list that records the levels in
##'         each factor variable.
##'     \item \code{contrasts}: Contrasts specified and used for each
##'         factor variable.
##'     \item \code{convergCode}: \code{code} returned by function
##'         \code{\link[stats]{nlm}}, which is an integer indicating why the
##'         optimization process terminated. \code{help(nlm)} for details.
##'     \item \code{logL}: Log likelihood of the fitted model.
##'     \item \code{fisher}: Observed Fisher information matrix.
##' }
##'
##'
##' @references
##'
##' Wang, W., Chen, K., & Yan, J. (2017+).  Extended Cox Model by ECM Algorithm
##' for Uncertain Survival Records Due to Imperfect Data Integration. (working
##' in progress)
##'
##'
##'
##' @examples
##' library(intsurv)
##' intCox(survi())



## implementation of ECM algorithm to Cox model
intCox <- function(formula, data, subset, na.action, contrasts = NULL,
                   start = list(), control = list(), ...) {

    ## record the function call to return
    Call <- match.call()

    ## arguments check
    if (missing(formula))
        stop("Argument 'formula' is required.")
    if (missing(data))
        data <- environment(formula)
    if (! with(data, inherits(eval(formula[[2L]]), "survi")))
        stop("Response in formula must be a 'surve' object.")

    ## Prepare data: ID, time, event ~ X(s)
    mcall <- match.call(expand.dots = FALSE)
    mmcall <- match(c("formula", "data", "subset", "na.action"),
                    names(mcall), 0L)
    mcall <- mcall[c(1L, mmcall)]

    ## drop unused levels in factors
    mcall$drop.unused.levels <- TRUE
    mcall[[1L]] <- quote(stats::model.frame)
    mf <- eval(mcall, parent.frame())
    mt <- attr(mf, "terms")
    mm <- stats::model.matrix(formula, mf, contrasts.arg = contrasts)

    ## number of covariates excluding intercept
    if ((nBeta <- ncol(mm) - 1L) <= 0)
        stop("Covariates must be specified in formula.")

    ## covariates' names
    covar_names <- colnames(mm)[- 1L]

    ## data
    dat <- as.data.frame(cbind(mf[, 1L], mm[, - 1L]))
    colnames(dat) <- c("ID", "time", "event", covar_names)
    nObs <- nrow(dat)
    nBeta <- ncol(dat) - 3L
    dat$eventInd <- dat$event == 1L

    ## start' values for 'nlm'
    startList <- c(start, list(nBeta_ = nBeta, dat_ = dat))
    start <- do.call("intCox_start", startList)

    ## 'control' for 'nlm'
    control <- c(control, list(censorRate0_ = start$censorRate0))
    control <- do.call("intCox_control", control)

    ## indicator for subjects having multiple record
    dat$dupIdx <- with(dat, ID %in% unique(ID[duplicated(ID)]))

    ## define some variables for ease of computing
    incDat <- dat[(orderInc <- with(dat, order(time, ID))), ]
    dupVec <- duplicated(incDat$time)
    tied <- any(dupVec)
    incDat$firstIdx <- ! dupVec
    h0Dat <- h_cDat <- data.frame(time = incDat$time[incDat$firstIdx])
    xMat <- as.matrix(incDat[, 4L : (3L + nBeta)])

    ## from different starting values of piVec
    logL_max0 <- - Inf
    for (oneStart in start$censorRate) {

        if (is.null(start$piVec)) {
            ## determined by censorRate
            incDat$piVec <- piVec <- initPi(oneStart, dat = incDat)
        } else {
            incDat$piVec <- piVec <- start$piVec[orderInc]
        }

        ## trace the log-likelihood for observed data
        logL <- rep(NA, control$iterlimEm)
        ## trace beta estimates from each iteration of ECM
        betaMat <- matrix(NA, nrow = control$iterlimEm + 1L, ncol = nBeta)
        betaMat[1L, ] <- start$beta
        tolPi <- sqrt(control$tolEm)

        for (iter in seq_len(control$iterlimEm)) {
            oneFit <- oneECMstep(betaHat = betaMat[iter, ], h0Dat = h0Dat,
                                 h_cDat = h_cDat, dat = incDat, xMat = xMat,
                                 tied = tied, control = control)
            ## log likehood
            logL[iter] <- oneFit$logL

            ## always update p_jk or not? maybe yes
            if (control$alwaysUpdatePi || (iter > 1 && tol < tolPi))
                incDat$piVec <- oneFit$piVec

            ## update beta estimates
            betaEst <- oneFit$betaEst
            betaMat[iter + 1L, ] <- betaEst$estimate
            tol <- sum((betaMat[iter + 1L, ] - betaMat[iter, ]) ^ 2) /
                sum((betaMat[iter + 1L, ] + betaMat[iter, ]) ^ 2)
            if (tol < control$tolEm) {
                betaHat <- betaEst$estimate
                break
            }
        }

        ## keep the one fit maximizing observed log likelihood
        logL <- stats::na.omit(logL)
        logL_max <- logL[length(logL)]
        if (logL_max > logL_max0) {
            logL_max0 <- logL_max
            logL0 <- logL
            oneFit0 <- oneFit
            betaMat0 <- betaMat
            piVec0 <- piVec
            censorRate0 <- oneStart
        }
    }

    ## clean-up NA's
    betaMat0 <- stats::na.omit(betaMat0)
    attr(logL0, "na.action") <- attr(logL0, "class") <-
        attr(betaMat0, "na.action") <- attr(betaMat0, "class") <- NULL
    ## prepare for outputs
    piEst <- oneFit0$piVec[(reOrderIdx <- order(orderInc))]
    start$piVec <- piVec0[reOrderIdx]
    start$censorRate0 <- censorRate0
    ## update results
    h0Dat$h0Vec <- oneFit0$h0Vec
    h_cDat$h_cVec <- oneFit0$h_cVec
    incDat$xExp <- oneFit0$xExp
    incDat$piVec <- oneFit0$piVec
    betaEst <- oneFit0$betaEst
    betaHat <- as.numeric(betaEst$estimate)
    ## numerical approximation of I_oc fisher information matrix
    ## I_oc <- approxIoc(dat = incDat, xMat = xMat, tied = tied, nIter = 100)
    I_oc <- betaEst$hessian

    se_prem <- NA
    if (! control$noSE) {
        ## dm matrix
        ## incDat$piVec <- initPi(censorRate0, dat = incDat, equally = FALSE)
        ## tmpRes <- oneECMstep(betaHat, h0Dat = h0Dat, h_cDat = h_cDat,
        ##                      dat = incDat, xMat = xMat, tied = tied,
        ##                      control = control)
        ## incDat$piVec <- tmpRes$piVec
        dmMat <- dmECM(betaEst = betaHat, h0Dat = h0Dat, h_cDat = h_cDat,
                       dat = incDat, xMat = xMat, tied = tied,
                       control = control)

        ## variance-covariance matrix by SECM
        invI_oc <- solve(I_oc)
        secmVar <- invI_oc + invI_oc %*% dmMat %*% solve(diag(1, nBeta) - dmMat)
        secmVar <- (secmVar + t(secmVar)) / 2
        ## secmVar <- solve((diag(1, nBeta) - dmMat) %*% I_oc)

        ## se estimatation by multiple imputation method
        ## miVar <- imputeVar(incDat, xMat, nImpute = 30)

        ## se estimatation by PRES, FIXME
        ## I_oMat <- I_o(betaEst = betaHat, h0Dat = h0Dat, h_cDat = h_cDat,
        ##               dat = incDat, xMat = xMat, tied = tied,
        ##               control = control)
        ## invI_o <- solve(I_oMat)
        se_prem <- as.numeric(sqrt(diag(secmVar)))
    }

    ## estimates for beta
    est_beta <- matrix(NA, nrow = nBeta, ncol = 5L)
    colnames(est_beta) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
    rownames(est_beta) <- covar_names
    est_beta[, 1L] <- betaHat
    est_beta[, 2L] <- exp(est_beta[, "coef"])
    ## est_beta[, 3L] <- NA
    ## est_beta[, 4L] <- est_beta[, 1L] / est_beta[, 3L]
    ## est_beta[, 5L] <- 2 * stats::pnorm(- abs(est_beta[, 5L]))

    ## output: na.action
    na.action <- if (is.null(attr(mf, "na.action")))
                     options("na.action")[[1L]]
                 else
                     paste("na.", class(attr(mf, "na.action")))

    ## output: contrasts
    contrasts <- if (is.null(contrasts))
                     list(contrasts = NULL)
                 else
                     attr(mm, "contrasts")

    ## results to return
    results <- methods::new("intCox",
                            call = Call,
                            formula = formula,
                            data = as.data.frame(data),
                            nObs = nObs,
                            estimates = list(beta = est_beta,
                                             pi = piEst,
                                             h0 = h0Dat),
                            control = control,
                            start = start,
                            na.action = na.action,
                            xlevels = .getXlevels(mt, mf),
                            contrasts = contrasts,
                            convergCode = betaEst$code,
                            partLogL = - betaEst$minimum,
                            logL = logL0,
                            fisher = betaEst$hessian)
    ## return
    results
}


### internal functions =========================================================
## one iteration for one row of DM matrix
oneRowDM <- function(ind, betaEst, h0Dat, h_cDat, dat, xMat, tied, control)
{
    h <- control$h
    nBeta <- length(betaEst)
    baseVec <- rep(0, nBeta)
    baseVec[ind] <- h
    theta1 <- betaEst - 2 * baseVec
    theta2 <- betaEst - baseVec
    theta3 <- betaEst + baseVec
    theta4 <- betaEst + 2 * baseVec
    oneFit1 <- oneECMstep(theta1, h0Dat, h_cDat, dat, xMat, tied, control)
    oneFit2 <- oneECMstep(theta2, h0Dat, h_cDat, dat, xMat, tied, control)
    oneFit3 <- oneECMstep(theta3, h0Dat, h_cDat, dat, xMat, tied, control)
    oneFit4 <- oneECMstep(theta4, h0Dat, h_cDat, dat, xMat, tied, control)
    beta1 <- oneFit1$betaEst$estimate
    beta2 <- oneFit2$betaEst$estimate
    beta3 <- oneFit3$betaEst$estimate
    beta4 <- oneFit4$betaEst$estimate
    (beta1 - beta4 + 8 * (beta3 - beta2)) / (12 * h)
}


## DM matrix for ECM by PREM
dmECM <- function(betaEst, h0Dat, h_cDat, dat, xMat, tied, control)
{
    dmList <- lapply(seq_along(betaEst), oneRowDM,
                     betaEst = betaEst, h0Dat = h0Dat, h_cDat = h_cDat,
                     dat = dat,  xMat = xMat, tied = tied, control = control)
    do.call(rbind, dmList)
}


## Fisher score function
scoreFun <- function(betaEst, h0Dat, h_cDat, dat, xMat, tied, control)
{
    ## update results involving beta estimates
    dat$betaX <- as.numeric(xMat %*% betaEst)
    dat$xExp <- exp(dat$betaX)
    dat$xExp <- ifelse(is.infinite(dat$xExp), 1e50, dat$xExp)

    ## use initial value of pi
    dat$p_jk <- dat$piVec
    h0Dat$h0Vec <- h0t(dat, tied = tied)
    h_cDat$h_cVec <- h_c(dat, tied = tied)

    ## update baseline hazard rate of event times
    h0Dat$H0Vec <- cumsum(h0Dat$h0Vec)
    ## update baseline hazard rate of censoring times
    h_cDat$H_cVec <- cumsum(h_cDat$h_cVec)

    if (tied) {
        time_idx <- match(dat$time, h0Dat$time)
        dat$h0Vec <- with(dat, ifelse(eventInd, h0Dat$h0Vec[time_idx], 0))
        dat$H0Vec <- h0Dat$H0Vec[time_idx]
        dat$h_cVec <- with(dat, ifelse(eventInd, 0, h_cDat$h_cVec[time_idx]))
        dat$H_cVec <- h_cDat$H_cVec[time_idx]
    } else {
        dat$h0Vec <- with(dat, ifelse(eventInd, h0Dat$h0Vec, 0))
        dat$H0Vec <- h0Dat$H0Vec
        dat$h_cVec <- with(dat, ifelse(eventInd, 0, h_cDat$h_cVec))
        dat$H_cVec <- h_cDat$H_cVec
    }

    dat$hVec <- with(dat, h0Vec * xExp)
    dat$HVec <- with(dat, H0Vec * xExp)

    ## scaling
    foo <- sum(dat$HVec)
    dat$sVec <- exp(- dat$HVec + log(foo)) / foo
    foo <- sum(dat$H_cVec)
    dat$G_cVec <- exp(- dat$H_cVec + log(foo)) / foo

    ## compute p_jk for each subject
    ## for observed log-likelihood function
    dat$p_jk_numer <- with(dat, ifelse(eventInd,
                                       piVec * hVec * sVec * G_cVec,
                                       piVec * sVec * h_cVec * G_cVec))

    s_p_jk_denom <- with(dat, tapply(p_jk_numer, ID, FUN = sum))
    idx <- match(as.character(dat$ID), names(s_p_jk_denom))
    dat$p_jk_denom <- s_p_jk_denom[idx]
    dat$p_jk <- with(dat, ifelse(dupIdx,
                          ifelse(p_jk_denom == 0, piVec,
                                 p_jk_numer / p_jk_denom), 1))

    ## building blocks
    parSeq <- seq_along(betaEst)
    xMatDeltaN <- xMat[dat$eventInd, ] * dat[dat$eventInd, "p_jk"]
    delta_tildeN <- deltaTildeN(dat, tied)
    k_0 <- k0(dat, tied)
    k_1 <- k1(parSeq, dat, xMat, tied)
    dLbeta(xMatDeltaN, k_0, k_1, delta_tildeN)
}


## SE from PRES
oneRowI_o <- function(ind, betaEst, h0Dat, h_cDat, dat, xMat, tied, control)
{
    h <- control$h
    nBeta <- length(betaEst)
    baseVec <- rep(0, nBeta)
    baseVec[ind] <- h
    theta1 <- betaEst - 2 * baseVec
    theta2 <- betaEst - baseVec
    theta3 <- betaEst + baseVec
    theta4 <- betaEst + 2 * baseVec
    oneFit1 <- scoreFun(theta1, h0Dat, h_cDat, dat, xMat, tied, control)
    oneFit2 <- scoreFun(theta2, h0Dat, h_cDat, dat, xMat, tied, control)
    oneFit3 <- scoreFun(theta3, h0Dat, h_cDat, dat, xMat, tied, control)
    oneFit4 <- scoreFun(theta4, h0Dat, h_cDat, dat, xMat, tied, control)
    (oneFit1 - oneFit4 + 8 * (oneFit3 - oneFit2)) / (12 * h)
}


## DM matrix for ECM by PREM
I_o <- function(betaEst, h0Dat, h_cDat, dat, xMat, tied, control)
{
    resList <- lapply(seq_along(betaEst), oneRowI_o,
                      betaEst = betaEst, h0Dat = h0Dat, h_cDat = h_cDat,
                      dat = dat,  xMat = xMat, tied = tied, control = control)
    do.call(rbind, resList)
}


## sample latent indicators based on estiamted posterior prob.
rLatent <- function(dat) {
    tmpList <- with(dat, tapply(piVec, ID, function(a) {
        if (sum(a))
            return(rmultinom(1L, size = 1L, prob = a))
        ## else return
        0
    }))
    unlist(tmpList)
}


## one imputation
oneImpute <- function(dat, xMat) {
    ## one latent sample
    dat$latent <- rLatent(dat)
    uniDat <- subset(dat, latent == 1L)
    oneFit <- survival::coxph(survival::Surv(time, event) ~ xMat,
                              data = dat, ties = "breslow")
    est <- oneFit$coefficients
    betaVar <- diag(oneFit$var)
    c(est, betaVar)
}


## se estimates by multiple imputation
imputeVar <- function(dat, xMat, nImpute) {
    nBeta <- ncol(xMat)
    idx <- seq_len(nBeta)
    coefMat <- replicate(nImpute, oneImpute(dat, xMat))
    eVar <- apply(coefMat[idx, ], 1, var)
    adj <- 1 + 1 / nImpute
    mVar <- rowMeans(coefMat[- idx, ])
    mVar + adj * eVar
}


## I_oc matrix from one latent sample
oneIoc <- function(parSeq, dat, xMat, tied) {
    ## one latent sample
    dat$p_jk <- rLatent(dat)
    ## prepare for fisher information matrix
    delta_tildeN <- deltaTildeN(dat, tied)
    k_0 <- k0(dat, tied)
    k_1 <- k1(parSeq, dat, xMat, tied)
    k_2 <- k2(parSeq, dat, xMat, tied)
    ## one fisher information matrix in vector
    - as.vector(d2Lbeta(parSeq, k_0, k_1, k_2, delta_tildeN))
}


## approximation of I_oc matrix
approxIoc <- function(dat, xMat, tied, nIter = 1e3) {
    nBeta <- ncol(xMat)
    parSeq <- seq_len(nBeta)
    tmpMat <- replicate(nIter, oneIoc(parSeq, dat, xMat, tied))
    tmpVec <- rowMeans(tmpMat)
    matrix(tmpVec, ncol = nBeta)
}


## perform one step of EM algorithm
oneECMstep <- function(betaHat, h0Dat, h_cDat, dat, xMat, tied, control)
{
    ## update results involving beta estimates
    dat$betaX <- as.numeric(xMat %*% betaHat)
    dat$xExp <- exp(dat$betaX)
    dat$xExp <- ifelse(is.infinite(dat$xExp), 1e50, dat$xExp)

    ## use initial value of pi
    dat$p_jk <- dat$piVec
    h0Dat$h0Vec <- h0t(dat, tied = tied)
    h_cDat$h_cVec <- h_c(dat, tied = tied)

    ## E-step ------------------------------------------------------------------
    ## update baseline hazard rate of event times
    h0Dat$H0Vec <- cumsum(h0Dat$h0Vec)
    ## update baseline hazard rate of censoring times
    h_cDat$H_cVec <- cumsum(h_cDat$h_cVec)

    if (tied) {
        time_idx <- match(dat$time, h0Dat$time)
        dat$h0Vec <- with(dat, ifelse(eventInd, h0Dat$h0Vec[time_idx], 0))
        dat$H0Vec <- h0Dat$H0Vec[time_idx]
        dat$h_cVec <- with(dat, ifelse(eventInd, 0, h_cDat$h_cVec[time_idx]))
        dat$H_cVec <- h_cDat$H_cVec[time_idx]
    } else {
        dat$h0Vec <- with(dat, ifelse(eventInd, h0Dat$h0Vec, 0))
        dat$H0Vec <- h0Dat$H0Vec
        dat$h_cVec <- with(dat, ifelse(eventInd, 0, h_cDat$h_cVec))
        dat$H_cVec <- h_cDat$H_cVec
    }

    dat$hVec <- with(dat, h0Vec * xExp)
    dat$HVec <- with(dat, H0Vec * xExp)
    dat$sVec <- exp(- dat$HVec)
    dat$G_cVec <- exp(- dat$H_cVec)

    ## compute p_jk for each subject
    ## for observed log-likelihood function
    dat$p_jk_numer <- with(dat, ifelse(eventInd,
                                       piVec * hVec * sVec * G_cVec,
                                       piVec * sVec * h_cVec * G_cVec))

    s_p_jk_denom <- with(dat, tapply(p_jk_numer, ID, FUN = sum))
    idx <- match(as.character(dat$ID), names(s_p_jk_denom))
    dat$p_jk_denom <- s_p_jk_denom[idx]
    dat$p_jk <- with(dat, ifelse(dupIdx,
                          ifelse(p_jk_denom == 0, piVec,
                                 p_jk_numer / p_jk_denom), 1))

    ## help speed up convergence?
    ## dat$p_jk[dat$p_jk < .Machine$double.eps] <- 0
    ## dat$p_jk[dat$p_jk > 1 - .Machine$double.eps] <- 1

    ## CM-steps ----------------------------------------------------------------
    ## update beta
    betaEst <- stats::nlm(logLbeta, p = betaHat, dat = dat, xMat = xMat,
                          tied = tied, hessian = TRUE, check.analyticals = TRUE,
                          gradtol = control$gradtol, stepmax = control$stepmax,
                          steptol = control$steptol, iterlim = control$iterlim)

    ## log-likelihood function under observed data
    logL <- sum(log(dat$p_jk_denom))

    ## update h0_jk and h_c_jk with previous (or initial) estimates of beta
    ## h0Vec <- h0t(dat, tied)
    ## h_cVec <- h_c(dat, tied)
    list(betaEst = betaEst, h0Vec = h0Dat$h0Vec, h_cVec = h_cDat$h_cVec,
         piVec = dat$p_jk, logL = logL, xExp = dat$xExp)
}


## profile log-likelihood function of beta
logLbeta <- function(param, dat, xMat, tied)
{
    ## update retults depends on beta estimates, param
    dat$betaX <- as.vector(xMat %*% param)
    dat$xExp <- exp(dat$betaX)
    dat$xExp <- ifelse(is.infinite(dat$xExp), 1e50, dat$xExp)

    ## prepare intermediate results for later computation
    parSeq <- seq_along(param)
    xMatDeltaN <- xMat[dat$eventInd, ] * dat[dat$eventInd, "p_jk"]
    delta_tildeN <- deltaTildeN(dat, tied)
    betaXdeltaN <- with(dat, eventInd * betaX * p_jk)

    k_0 <- k0(dat, tied)
    k_1 <- k1(parSeq, dat, xMat, tied)
    k_2 <- k2(parSeq, dat, xMat, tied)

    ## profile log-likelihood of beta in EM
    pell <- sum(betaXdeltaN) - sum(ifelse(delta_tildeN > 0,
                                          delta_tildeN * log(k_0), 0))

    ## penalty term to avoid solution with xExp being Inf
    penal_inf <- any(is.infinite(dat$xExp)) * 1e20
    negLogL <- - pell + penal_inf

    ## gradient
    gradLogL <- dLbeta(xMatDeltaN, k_0, k_1, delta_tildeN)
    attr(negLogL, "gradient") <- - gradLogL
    ## hessian
    hesMat <- d2Lbeta(parSeq, k_0, k_1, k_2, delta_tildeN)
    attr(negLogL, "hessian") <- - hesMat
    ## return
    negLogL
}


## compute baseline hazard rate function
h0t <- function(dat, tied) {
    numer <- deltaTildeN(dat, tied)
    denom <- k0(dat, tied)
    ifelse(numer > 0, numer / denom, 0)
}


## building blocks that follows notation in manuscript
deltaTildeN <- function(dat, tied) {
    out <- with(dat, eventInd * p_jk)
    if (tied)
        out <- tapply(out, dat$time, FUN = sum)
    as.numeric(out)
}


## baseline hazard rate for censoring time
h_c <- function(dat, tied) {
    numer <- deltaC(dat, tied)
    denom <- rev(cumsum(rev(dat$p_jk)))
    if (tied)
        denom <- denom[dat$firstIdx]
    ifelse(numer > 0, numer / denom, 0)
}


deltaC <- function(dat, tied) {
    out <- with(dat, (! eventInd) * p_jk)
    if (tied)
        out <- tapply(out, dat$time, FUN = sum)
    as.numeric(out)
}


k0 <- function(dat, tied) {
    p_jk_xExp <- with(dat, p_jk * xExp)
    out <- rev(cumsum(rev(p_jk_xExp)))
    if (tied)
        out <- out[dat$firstIdx]
    out
}


k1 <- function(parSeq, dat, xMat, tied) {
    ## matrix of dimension # unique event time by # parameters
    out <- sapply(parSeq, function(ind, dat, xMat) {
        p_jk_xExp_x <- with(dat, p_jk * xExp) * xMat[, ind]
        rev(cumsum(rev(p_jk_xExp_x)))
    }, dat = dat, xMat = xMat)
    if (tied)
        out <- out[dat$firstIdx, ]
    out
}


k2 <- function(parSeq, dat, xMat, tied) {
    ind_grid <- expand.grid(parSeq, parSeq)
    ## matrix of dimension # unique event time by (# parameters) ^ 2
    out <- mapply(function(ind1, ind2) {
        p_jk_xExp_x2 <- with(dat, p_jk * xExp) * xMat[, ind1] * xMat[, ind2]
        rev(cumsum(rev(p_jk_xExp_x2)))
    }, ind_grid[, 1L], ind_grid[, 2L])
    if (tied)
        out <- out[dat$firstIdx, ]
    out
}


dLbeta <- function(xMatDeltaN, k_0, k_1, delta_tildeN) {
    sum_jk <- colSums(xMatDeltaN)
    int_t <- colSums(na.omit(k_1 / k_0 * delta_tildeN))
    sum_jk - int_t
}


d2Lbeta <- function(parSeq, k_0, k_1, k_2, delta_tildeN) {
    ## part 1
    nPar <- length(parSeq)
    mat1 <- na.omit(k_2 / k_0 * delta_tildeN)
    part1 <- matrix(colSums(mat1), nPar, nPar)

    ## part 2
    mat2 <- k_1 / k_0
    mat2[is.na(mat2)] <- 0
    parGrid <- expand.grid(parSeq, parSeq)
    mat2plus <- mapply(function(ind1, ind2) {
        mat2[, ind1] * mat2[, ind2]
    }, parGrid[, 1L], parGrid[, 2L])
    part2 <- matrix(colSums(mat2plus * delta_tildeN), nPar, nPar)
    ## return
    part2 - part1
}


## determine initial piVec from given censorRate
initPi <- function(censorRate, dat, equally = FALSE, ...) {
    ## mixture probability for each subject: piVec
    numTab <- table(dat$ID)
    dat$numRecord <- numTab[match(as.character(dat$ID), names(numTab))]

    ## subject ID with censoring records
    cenID <- with(subset(dat, ! eventInd), unique(ID))
    cenIdx <- as.integer(dat$ID %in% cenID)

    ## for subjects with multiple records
    idx1 <- with(dat, numRecord > 1L & eventInd)
    idx2 <- with(dat, numRecord > 1L & (! eventInd))

    piVec <- rep(1L, NROW(dat))
    if (equally) {
        piVec[idx1 | idx2] <- 1 / dat[idx1 | idx2, "numRecord"]
    } else {
        piVec[idx1] <- (1 - censorRate * cenIdx[idx1]) /
            (dat[idx1, "numRecord"] - cenIdx[idx1])
        piVec[idx2] <- censorRate
    }
    piVec
}


intCox_start <- function(beta, censorRate, piVec, ..., nBeta_, dat_)
{
    dupID <- with(dat_, unique(ID[duplicated(ID)]))
    uniDat <- base::subset(dat_, ! ID %in% dupID)
    censorRate0 <- round(1 - mean(uniDat$event), 2)
    if (missing(censorRate)) {
        ## set censorRate from sample truth data
        ## if missing at random, the true censoring rate
        ## can be estimated by true data of unique records
        step_by <- 0.02
        censorRate <- seq.int(0, 1, step_by)
    } else if (any(censorRate > 1 | censorRate < 0))
        stop(paste("Starting prob. of censoring case being true",
                   "should between 0 and 1."))

    ## initialize covariate coefficient: beta
    if (missing(beta)) {
        ## if high censoring for subjects having unique records
        if (mean(uniDat$event) < 0.01) {
            beta <- rep(0, nBeta_)
        } else {
            uniDat$eventInd <- NULL
            tmp <- tryCatch(
                survival::coxph(survival::Surv(time, event) ~
                                    as.matrix(uniDat[, - seq_len(3L)]),
                                data = uniDat),
                warning = function(w) {
                    warning(w)
                    return(NULL)
                })
            beta <- if (is.null(tmp))
                        rep(0, nBeta_)
                    else
                        as.numeric(tmp$coefficients)
        }
    } else {
        beta <- as.numeric(beta)
        if (length(beta) != nBeta_)
            stop(paste("Number of starting values for coefficients of",
                       "covariates does not match with the specified formula."))
    }

    if (missing(piVec)) {
        piVec <- NULL
    } else {
        if (length(piVec) != nrow(dat_))
            stop("'piVec' must have same length with number of rows of data.")
        if (any(piVec > 1 | piVec < 0))
            stop("'piVec' has to be between 0 and 1.")
        censorRate <- NA                # will not be used
    }

    ## return
    list(beta = beta, censorRate = censorRate,
         piVec = piVec, censorRate0 = censorRate0)
}


intCox_control <- function(gradtol = 1e-6, stepmax = 1e2,
                         steptol = 1e-6, iterlim = 1e2,
                         tolEm = 1e-4, iterlimEm = 1e2,
                         h = sqrt(tolEm), alwaysUpdatePi = NULL, ...,
                         censorRate0_, noSE_ = TRUE)
{
    ## controls for function stats::nlm
    if (! is.numeric(gradtol) || gradtol <= 0)
        stop("value of 'gradtol' must be > 0.")
    if (! is.numeric(stepmax) || stepmax <= 0)
        stop("value of 'stepmax' must be > 0.")
    if (! is.numeric(steptol) || steptol <= 0)
        stop("value of 'steptol' must be > 0.")
    if (! is.numeric(iterlim) || iterlim <= 0)
        stop("maximum number of iterations must be > 0.")

    ## determining convergence of EM
    if (! is.numeric(tolEm) || tolEm <= 0)
        stop("value of 'tolEm' must be > 0.")
    if (! is.numeric(iterlimEm) || iterlimEm <= 0)
        stop("maximum number of iterations for EM must be > 0.")

    ## determining the computation of DM matrix
    if (! is.numeric(h) || h < 0)
        stop("'h' has to be a positive number.")

    ## automatically determine whether always update pi's
    if (is.null(alwaysUpdatePi))
        alwaysUpdatePi <- ifelse(censorRate0_ < 0.8, TRUE, FALSE)

    ## return
    list(gradtol = gradtol, stepmax = stepmax,
         steptol = steptol, iterlim = iterlim,
         tolEm = tolEm, iterlimEm = iterlimEm,
         h = h, alwaysUpdatePi = alwaysUpdatePi,
         noSE_ = noSE_)
}
