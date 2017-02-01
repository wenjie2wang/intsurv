################################################################################
### functions fitting Cox model with uncertain survival record data
### version controlled by git
################################################################################


##' FIXME
##'
##'
##' @examples
##' ## temp code for testing only
##' library(survival)
##' source("class.R")
##' source("coef.R")
##' source("fit.R")
##' source("../simulation/simuData.R")
##' source("../simulation/simuFun.R")
##'
##' set.seed(1216)
##' dat <- simuWeibull(nSubject = 1000,
##'                    maxNum = 1, nRecordProb = c(0.7, 0.3),
##'                    matchCensor = 4 / 15, matchEvent = 4 / 5,
##'                    censorMax = 12.5, censorMin = 0.5,
##'                    lambda = 0.045, rho = 0.7,
##'                    fakeLambda1 = 0.045 * exp(- 3),
##'                    fakeLambda2 = 0.045 * exp(3),
##'                    mixture = 0.5, eventOnly = FALSE)
##' ## dat$obsTime <- round(dat$obsTime, 2)
##' temp <- coxEm(Surve(ID, obsTime, eventInd) ~ x1 + x2,
##'               data = dat, control = list(alwaysUpdatePi = FALSE),
##'               start = list(censorRate = seq.int(0, 1, 0.05)))
##' temp@logL
##' temp@start$censorRate0
##' summar(list(temp))
##'
##' temp <- coxEm(Surve(ID, obsTime, eventInd) ~ x1 + x2, data = dat,
##'               control = list(alwaysUpdatePi = TRUE),
##'               start = list(censorRate = seq.int(0, 0.95, 0.01)))
##' temp@logL
##' temp@start$censorRate0
##' summar(list(temp))
##'
##' piVec0 <- with(dat, as.integer(! duplicated(ID)))
##' tmp <- coxEm(Surve(ID, obsTime, eventInd) ~ x1 + x2, data = dat,
##'              start = list(beta = c(0, 0), piVec = piVec0, censorRate = 1))
##' summar(list(tmp))
##' tmp@logL
##' ## temp <- coxEm(Surve(ID, obsTime, eventInd) ~ x1 + x2, data = dat,
##' ##               start = list(beta = c(1, 1)))
##'
##' ## temp@estimates$beta
##' tmpDat <- cbind(dat, piEst = round(temp@estimates$piEst, 3))
##' dupID <- with(tmpDat, unique(ID[duplicated(ID)]))
##' subset(tmpDat, ID %in% dupID)
##' xtabs(~ eventInd + piEst, tmpDat, latentInd != 1L)
##' xtabs(~ eventInd + piEst, tmpDat, latentInd == 1L)
##'
##' summar(list(temp))
##' naiveCox(temp)
##' uniOnlyCox(temp)
##' oracleCox(temp)
##' oracleWb(temp, rho0 = 0.7)
##'
##' ## test on the true data of unique records
##' trueDat <- dat[with(dat, ! duplicated(ID)), ]
##' tmp <- coxEm(Surve(ID, obsTime, eventInd) ~ x1 + x2, data = trueDat,
##'              start = list(beta = c(0, 0)))
##' tmp@estimates$beta
##'
##' ## test on tied event times
##' test1 <- list(ID = seq_len(7),
##'               time = c(4, 3, 1, 1, 2, 2, 3),
##'               status = c(1, 1, 1, 0, 1, 1, 0),
##'               x = c(0, 2, 1, 1, 1, 0, 0),
##'               sex = c(0, 0, 0, 0, 1, 1, 1))
##' coxph(Surv(time, status) ~ x + sex, test1, ties = "breslow")
##' coef(coxEm(Surve(ID, time, status) ~ x + sex, test1,
##'            start = list(beta = c(0, 0))))
##'
##' ## @import data.table? for faster aggregation?


## implementation of ECM algorithm to Cox model
coxEm <- function(formula, data, subset, na.action, contrasts = NULL,
                  start = list(), control = list(), ...) {

    ## record the function call to return
    Call <- match.call()

    ## arguments check
    if (missing(formula))
        stop("Argument 'formula' is required.")
    if (missing(data))
        data <- environment(formula)
    if (! with(data, inherits(eval(Call[[2L]][[2L]]), "Surve")))
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
    mm <- stats::model.matrix(formula, data = mf, contrasts.arg = contrasts)

    ## number of covariates excluding intercept
    if ((nBeta <- ncol(mm) - 1L) <= 0)
        stop("Covariates must be specified in formula.")

    ## covariates' names
    covar_names <- colnames(mm)[- 1L]

    ## data
    dat <- as.data.frame(cbind(mf[, 1L][, seq_len(3L)], mm[, - 1L]))
    colnames(dat) <- c("ID", "time", "event", covar_names)
    nObs <- nrow(dat)
    nBeta <- ncol(dat) - 3L
    dat$eventInd <- dat$event == 1L

    ## 'control' for 'nlm'
    control <- do.call("coxEmControl", control)

    ## start' values for 'nlm'
    startList <- c(start, list(nBeta_ = nBeta, dat_ = dat))
    start <- do.call("coxEmStart", startList)

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
                                 tied = tied,
                                 control = control, start = start)
            ## log likehood
            logL[iter] <- oneFit$logL

            ## always update p_jk or not? maybe yes
            if (control$alwaysUpdatePi || (iter > 1 && tol < tolPi))
                incDat$piVec <- oneFit$piVec

            ## update beta estimates
            betaEst <- oneFit$betaEst
            betaMat[iter + 1L, ] <- betaEst$estimate
            tol <- max(abs((betaMat[iter + 1L, ] - betaMat[iter, ]) /
                           (betaMat[iter + 1L, ] + betaMat[iter, ])))
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
    ## update baseline hazard rate h0
    h0Dat$h0Vec <- oneFit0$h0Vec
    ## update baseline hazard rate of censoring times
    h_cDat$h_cVec <- oneFit0$h_cVec

    ## ## numerical approximation of I_oc fisher information matrix
    ## incDat$xExp <- oneFit$xExp
    ## I_oc <- approxIoc(dat = incDat, xMat = xMat, tied = tied, nIter = 100)

    ## ## DM matrix
    ## dmMat <- dmECM(betaMat = betaMat, betaEst = betaHat, h0Dat = h0Dat,
    ##                h_cDat = h_cDat, dat = incDat, xMat = xMat, tied = tied,
    ##                control = control)

    ## ## variance-covariance matrix by SECM
    ## invI_oc <- solve(I_oc)
    ## secmVar <- invI_oc + invI_oc %*% dmMat %*% solve(diag(1, nBeta) - dmMat)
    ## secmVar <- solve((diag(1, nBeta) - dmMat) %*% I_oc)

    ## se estimatation by multiple imputation method
    ## miVar <- imputeVar(incDat, xMat, nImpute = 30)

    ## estimates for beta
    ## est_beta <- matrix(NA, nrow = nBeta, ncol = 6L)
    est_beta <- matrix(NA, nrow = nBeta, ncol = 3L)
    ## colnames(est_beta) <- c("coef", "exp(coef)", "se(coef)",
    ##                         "se_SEM", "z", "Pr(>|z|)")
    ## colnames(est_beta) <- c("coef", "se_comp", "se_SEM", "se_MI")
    colnames(est_beta) <- c("coef", "se_comp", "se_boot")
    rownames(est_beta) <- covar_names
    betaEst <- oneFit0$betaEst
    se_vec <- sqrt(diag(solve(betaEst$hessian)))
    est_beta[, 1L] <- as.vector(betaEst$estimate)
    ## est_beta[, 2L] <- exp(est_beta[, 1L])
    est_beta[, 2L] <- as.vector(se_vec)
    est_beta[, 3L] <- NA
    ## est_beta[, 3L] <- as.vector(sqrt(diag(secmVar)))
    ## est_beta[, 4L] <- as.vector(sqrt(miVar))
    ## est_beta[, 6L] <- est_beta[, 1L] / est_beta[, 3L]
    ## est_beta[, 7L] <- 2 * stats::pnorm(- abs(est_beta[, 5L]))

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

    ## output: data
    colnames(dat) <- c(as.character(formula[[2]])[- 1L], covar_names)

    ## results to return
    results <- methods::new("coxEm",
                            call = Call,
                            formula = formula,
                            data = dat,
                            nObs = nObs,
                            estimates = list(beta = est_beta,
                                             piEst = piEst,
                                             h0Dat = h0Dat),
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
oneRowDM <- function(ind, betaOld, betaEst, h0Dat, h_cDat,
                     dat, xMat, tied, control) {
    betaOld_ind <- betaEst
    betaOld_ind[ind] <- betaOld[ind]
    oneFit <- oneECMstep(betaOld_ind, h0Dat, h_cDat, dat, xMat, tied, control)
    betaNew_ind <- oneFit$betaEst$estimate
    denom <- betaOld[ind] - betaEst[ind]
    numer <- betaNew_ind - betaEst
    numer / denom
}


## DM matrix for ECM
dmECM <- function(betaMat, betaEst, h0Dat, h_cDat, dat, xMat, tied, control) {
    ## set up iteration
    nBeta <- length(betaEst)
    idx <- seq_len(nBeta)
    transDM <- transDM_old <- matrix(0, nBeta, nBeta)
    tolVec <- rep(1, nBeta)
    iterMax <- nrow(betaMat)

    for (iter in min(iterMax, seq_len(control$iterlimSem))) {
        ## transpose of DM matrix
        transDM[, idx] <- sapply(idx, oneRowDM, betaOld = betaMat[iter, ],
                                 betaEst = betaEst, h0Dat = h0Dat,
                                 h_cDat = h_cDat, dat = dat, xMat = xMat,
                                 tied = tied, control = control)

        ## NaN produced when convergence of ECM
        if (any(is.nan(transDM))) {
            warning("DM matrix does not converge under given tolerance.")
            return(t(transDM_old))
        }

        ## determine convergence on each row of DM matrix (column of DM^T)
        dmArray <- array(c(transDM_old, transDM), dim = c(nBeta, nBeta, 2L))
        dmArray <- dmArray[idx, idx, seq_len(2), drop = FALSE]
        tolVec[idx] <- apply(dmArray, 2, function(vec2) {
            oneDM_old <- vec2[, 1L]
            oneDM_new <- vec2[, 2L]
            vecDiff <- oneDM_new - oneDM_old
            ## vecSum <- oneDM_new + oneDM_old
            ## max(abs(vecDiff / vecSum))
            max(abs(vecDiff))
        })
        if (iter > 1) {
            if (all(tolVec < control$tolSem)) break
            ## update for next iteration
            idx <- idx[tolVec > control$tolSem]
        }
        transDM_old <- transDM
    }
    t(transDM)
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

    ## help convergence?
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
initPi <- function(censorRate, dat, ...) {
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
    piVec[idx1] <- (1 - censorRate * cenIdx[idx1]) /
        (dat[idx1, "numRecord"] - cenIdx[idx1])
    piVec[idx2] <- censorRate
    ## piVec[idx1 | idx2] <- 1 / dat_[idx1 | idx2, "numRecord"]
    piVec
}


coxEmStart <- function(beta, censorRate, piVec, ..., nBeta_, dat_)
{
    dupID <- with(dat_, unique(ID[duplicated(ID)]))
    uniDat <- base::subset(dat_, ! ID %in% dupID)
    if (missing(censorRate)) {
        ## set censorRate from sample truth data
        ## if missing at random, the true censoring rate
        ## can be estimated by true data of unique records
        censorRate <- unique(c(seq.int(0, 1, 0.05),
                               1 - mean(uniDat$event)))
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
            tmp <- survival::coxph(survival::Surv(time, event) ~
                                       as.matrix(uniDat[, - seq_len(3L)]),
                                   data = uniDat)
            beta <- as.numeric(tmp$coefficients)
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
    list(beta = beta, censorRate = censorRate, piVec = piVec)
}


coxEmControl <- function(gradtol = 1e-6, stepmax = 1e2,
                         steptol = 1e-6, iterlim = 1e2,
                         tolEm = 1e-8, iterlimEm = 1e3,
                         tolSem = 1e-4, iterlimSem = 1e2,
                         alwaysUpdatePi = FALSE, ...)
{
    ## controls for function stats::nlm
    if (!is.numeric(gradtol) || gradtol <= 0)
        stop("value of 'gradtol' must be > 0")
    if (!is.numeric(stepmax) || stepmax <= 0)
        stop("value of 'stepmax' must be > 0")
    if (!is.numeric(steptol) || steptol <= 0)
        stop("value of 'steptol' must be > 0")
    if (!is.numeric(iterlim) || iterlim <= 0)
        stop("maximum number of iterations must be > 0")

    ## determining convergence of EM
    if (!is.numeric(tolEm) || tolEm <= 0)
        stop("value of 'tolEm' must be > 0")
    if (!is.numeric(iterlimEm) || iterlimEm <= 0)
        stop("maximum number of iterations for EM must be > 0")

    ## determining convergence of DM matrix in SEM
    if (!is.numeric(tolSem) || tolSem <= 0)
        stop("value of 'tolSem' must be > 0")
    if (!is.numeric(iterlimSem) || iterlimSem <= 0)
        stop("maximum number of iterations for SEM must be > 0")

    ## return
    list(gradtol = gradtol, stepmax = stepmax,
         steptol = steptol, iterlim = iterlim,
         tolEm = tolEm, iterlimEm = iterlimEm,
         tolSem = tolSem, iterlimSem = iterlimSem,
         alwaysUpdatePi = alwaysUpdatePi)
}
