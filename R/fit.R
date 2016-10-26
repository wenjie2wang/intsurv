################################################################################
### functions fitting Cox model with uncertain survival record data
### version controlled by git
################################################################################


##' FIXME
##'
##'
##' @examples
##' ## temp code for testing only
##'
##' source("class.R")
##' source("../simulation/simuData.R")
##'
##' set.seed(1216)
##' dat <- simuWeibull(nSubject = 200, maxNum = 2, nRecordProb = c(0.8, 0.2),
##'                    censorTime = 10, lambda = 1e-2, rho = 2, mixture = 0.5)
##' dat$obsTime <- round(dat$obsTime, 2)
##' temp <- coxEm(Surve(ID, obsTime, eventInd) ~ x1 + x2, data = dat)

## implementation of EM algorithm to Cox model
coxEm <- function(formula, data, subset, na.action, contrasts = NULL,
                  start = list(), control = list(), ...) {

    ## record the function call to return
    Call <- match.call()

    ## arguments check
    if (missing(formula))
        stop("Argument 'formula' is required.")
    if (missing(data))
        data <- environment(formula)
    if (! with(data, inherits(eval(Call[[2]][[2]]), "Surve")))
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
    dat <- as.data.frame(cbind(mf[, 1L][, seq_len(3)], mm[, - 1L]))
    colnames(dat) <- c("ID", "time", "event", covar_names)
    nObs <- nrow(dat)
    nBeta <- ncol(dat) - 3L

    ## define some variables for ease of computing
    incDat <- dat[(orderInc <- with(dat, order(time, ID))), ]
    xMat <- as.matrix(incDat[, - seq_len(3)])
    incDat$eventInd <- incDat$event == 1L

    ## 'control' for 'nlm'
    control <- do.call("coxEmControl", control)

    ## start' values for 'nlm'
    startList <- c(start, list(nBeta_ = nBeta, dat_ = incDat))
    start <- do.call("coxEmStart", startList)
    betaHat <- start$beta
    incDat$piVec <- start$piVec

    ## take care of possible ties
    h0Dat <- aggregate(start$h0Vec, list(time = incDat$time), FUN = sum)
    colnames(h0Dat)[2L] <- "h0Vec"

    logL <- rep(NA, control$iterlimEm)
    for (iter in seq_len(control$iterlimEm)) {
        oneFit <- oneEMstep(betaHat = betaHat, h0Dat = h0Dat,
                            dat = incDat, xMat = xMat, control = control)
        incDat$piVec <- oneFit$piVec
        h0Dat$h0Vec <- oneFit$h0Vec
        betaEst <- oneFit$betaEst
        betaHatNew <- betaEst$estimate
        ## log likehood
        logL[iter] <- oneFit$logL
        iter <- iter + 1L
        tol <- sqrt(sum((betaHatNew - betaHat) ^ 2)) / sqrt(sum(betaHat ^ 2))
        ## update beta estimates
        betaHat <- betaHatNew
        if (tol < control$tolEm) break
    }

    ## observed log likelihood from each iteration
    logL <- stats::na.omit(logL)
    attr(logL, "na.action") <- attr(logL, "class") <- NULL

    ## numerical approximation of I_oc fisher information matrix
    incDat$xExp <- oneFit$xExp
    I_oc <- approxIoc(dat = incDat, xMat = xMat, nIter = 100)

    ## DM matrix
    dmMat <- dm(betaOld = start$beta, betaEst = betaHat, h0Dat = h0Dat,
                dat = incDat, xMat = xMat, control = control)

    ## variance-covariance matrix by SEM
    invI_oc <- solve(I_oc)
    semVar <- invI_oc + invI_oc %*% dmMat %*% solve(diag(1, nBeta) - dmMat)

    ## estimates for beta
    est_beta <- matrix(NA, nrow = nBeta, ncol = 6L)
    colnames(est_beta) <- c("coef", "exp(coef)", "se(coef)",
                            "se_sem", "z", "Pr(>|z|)")
    rownames(est_beta) <- covar_names
    se_vec <- sqrt(diag(solve(betaEst$hessian)))
    est_beta[, 1] <- betaEst$estimate[seq_len(nBeta)]
    est_beta[, 2] <- exp(est_beta[, 1L])
    est_beta[, 3] <- se_vec[seq_len(nBeta)]
    est_beta[, 4] <- as.vector(sqrt(diag(semVar)))
    est_beta[, 5] <- est_beta[, 1L] / est_beta[, 3L]
    est_beta[, 6] <- 2 * stats::pnorm(- abs(est_beta[, 4L]))

    ## output: na.action
    na.action <- if (is.null(attr(mf, "na.action"))) {
                     options("na.action")[[1]]
                 } else {
                     paste("na.", class(attr(mf, "na.action")))
                 }
    ## output: contrasts
    contrasts <- if (is.null(contrasts)) {
                     list(contrasts = NULL)
                 } else {
                     attr(mm, "contrasts")
                 }
    ## results to return
    results <- methods::new("coxEm",
                            call = Call,
                            formula = formula,
                            data = dat,
                            nObs = nObs,
                            estimates = list(beta = est_beta,
                                             piVec = incDat$piVec,
                                             h0Dat = h0Dat,
                                             semVar = semVar),
                            control = control,
                            start = start,
                            na.action = na.action,
                            xlevels = .getXlevels(mt, mf),
                            contrasts = contrasts,
                            convergCode = betaEst$code,
                            partLogL = - betaEst$minimum,
                            logL = logL,
                            fisher = betaEst$hessian)
    ## return
    results
}


### internal functions =========================================================
## one iteration for one row of DM matrix
oneRowDM <- function(ind, betaOld, betaEst, h0Dat, dat, xMat, control){
    betaOld_ind <- betaEst
    betaOld_ind[ind] <- betaOld[ind]
    oneFit <- oneEMstep(betaOld_ind, h0Dat, dat, xMat, control)
    betaNew_ind <- oneFit$betaEst$estimate
    denom <- betaOld[ind] - betaEst[ind]
    numer <- betaNew_ind - betaEst
    numer / denom
}

## DM matrix for SEM
dm <- function(betaOld, betaEst, h0Dat, dat, xMat, control) {
    ## set up iteration
    nBeta <- length(betaEst)
    idx <- seq_len(nBeta)
    transDM <- transDM_old <- matrix(0, nBeta, nBeta)
    tolVec <- rep(1, nBeta)
    browser()
    for (iter in seq_len(control$iterlimSem)) {
        ## betaOld: \theta^{(t)}
        ## betaNew: \theta^{(t + 1)}
        oneFit <- oneEMstep(betaOld, h0Dat, dat, xMat, control)

        ## transpose of DM matrix
        transDM[, idx] <- sapply(idx, oneRowDM, betaOld = betaOld,
                                 betaEst = betaEst, h0Dat = h0Dat,
                                 dat = dat, xMat = xMat, control = control)

        ## determine convergence on each row of DM matrix (column of DM^T)
        dmArray <- array(c(transDM_old, transDM), dim = c(nBeta, nBeta, 2L))
        dmArray <- dmArray[idx, idx, seq_len(2), drop = FALSE]
        tolVec[idx] <- apply(dmArray, 2, function(vec2) {
            oneDM_old <- vec2[, 1L]
            oneDM_new <- vec2[, 2L]
            vecDiff <- oneDM_new - oneDM_old
            sum(vecDiff ^ 2) / sum(oneDM_old ^ 2)
        })
        if (all(tolVec < control$tolSem)) break

        ## update for next iteration
        betaOld <- oneFit$betaEst$estimate
        idx <- idx[tolVec > control$tolSem]
        transDM_old <- transDM
        dat$piVec <- oneFit$piVec
        h0Dat$h0Vec <- oneFit$h0Vec
    }

    t(transDM)
}


## sample latent indicators based on estiamted posterior prob.
rLatent <- function(dat) {
    tmpList <- aggregate(piVec ~ ID, data = dat, function(a) {
        rmultinom(1L, size = 1L, prob = a)
    })
    unlist(tmpList$piVec)
}

## I_oc matrix from one latent sample
oneIoc <- function(parSeq, dat, xMat) {
    ## one latent sample
    dat$p_jk <- rLatent(dat)
    ## prepare for fisher information matrix
    delta_tildeN <- deltaTildeN(dat)
    k_0 <- k0(dat)
    k_1 <- k1(parSeq, dat, xMat)
    k_2 <- k2(parSeq, dat, xMat)
    ## one fisher information matrix in vector
    - as.vector(d2Lbeta(parSeq, k_0, k_1, k_2, delta_tildeN))
}

## approximation of I_oc matrix
approxIoc <- function(dat, xMat, nIter = 1e3) {
    nBeta <- ncol(xMat)
    parSeq <- seq_len(nBeta)
    tmpMat <- replicate(nIter, oneIoc(parSeq, dat, xMat))
    tmpVec <- rowMeans(tmpMat)
    matrix(tmpVec, ncol = nBeta)
}

## perform one step of EM algorithm
oneEMstep <- function(betaHat, h0Dat, dat, xMat, control) {

    ## update results involving beta estimates
    dat$betaX <- as.vector(xMat %*% betaHat)
    dat$xExp <- pmax(exp(dat$betaX), .Machine$double.eps)

    ## update or initialize p_jk with piVec
    dat$p_jk <- dat$piVec

    ## update beta
    betaEst <- stats::nlm(logLbeta, p = betaHat, dat = dat, xMat = xMat,
                          hessian = TRUE, check.analyticals = TRUE,
                          gradtol = control$gradtol, stepmax = control$stepmax,
                          steptol = control$steptol, iterlim = control$iterlim)

    ## compute with updated beta estimates
    dat$betaX <- as.vector(xMat %*% betaEst$estimate)
    dat$xExp <- pmax(exp(dat$betaX), .Machine$double.eps)

    ## update baseline hazard rate
    h0Dat$H0Vec <- cumsum(h0Dat$h0Vec)
    time_idx <- match(dat$time, h0Dat$time)
    dat$h0Vec <- h0Dat$h0Vec[time_idx]
    dat$h0Vec <- with(dat, ifelse(eventInd, h0Vec, 0))
    dat$H0Vec <- h0Dat$H0Vec[time_idx]

    dat$hVec <- with(dat, h0Vec * xExp)
    dat$HVec <- with(dat, H0Vec * xExp)
    dat$sVec <- pmax(exp(- dat$HVec), .Machine$double.eps)

    ## compute p_jk for each subject
    dat$p_jk_numer <- with(dat, ifelse(eventInd, piVec * hVec * sVec,
                                       piVec * sVec))
    dat$p_jk_numer <- pmax(dat$p_jk_numer, .Machine$double.eps)
    p_jk_denom_dat <- aggregate(p_jk_numer ~ ID, data = dat, FUN = sum)
    idx <- match(dat$ID, p_jk_denom_dat$ID)
    dat$p_jk_denom <- p_jk_denom_dat[idx, "p_jk_numer"]
    dat$p_jk <- with(dat, p_jk_numer / p_jk_denom)

    ## log-likelihood function under observed data
    logL <- with(dat, sum(log(p_jk_denom[! duplicated(ID)])))

    ## update h0_jk with previous (or initial) estimates of beta
    h0Vec <- h0t(dat)
    list(betaEst = betaEst, h0Vec = h0Vec, piVec = dat$p_jk,
         logL = logL, xExp = dat$xExp)
}


## profile log-likelihood function of beta
logLbeta <- function(param, dat, xMat) {

    ## update retults depends on beta estimates, param
    dat$betaX <- as.vector(xMat %*% param)
    dat$xExp <- pmax(exp(dat$betaX), .Machine$double.eps)

    ## prepare intermediate results for later computation
    parSeq <- seq_along(param)
    xMatDeltaN <- xMat[dat$eventInd, ] * dat[dat$eventInd, "p_jk"]
    betaXdeltaN <- with(subset(dat, eventInd), betaX * p_jk)
    delta_tildeN <- deltaTildeN(dat)
    k_0 <- k0(dat)
    k_1 <- k1(parSeq, dat, xMat)
    k_2 <- k2(parSeq, dat, xMat)

    ## profile log-likelihood of beta in EM
    pell <- sum(betaXdeltaN) - sum(log(k_0) * delta_tildeN)
    ## penalty term to avoid solution with xExp being Inf
    penal <- any(dat$xExp == Inf) * 1e20
    negLogL <- - pell + penal

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
h0t <- function(dat) {
    numer <- deltaTildeN(dat)
    denom <- k0(dat)
    numer / denom
}

## building blocks that follows notation in manuscript
deltaTildeN <- function(dat) {
    tildeN_jk <- with(dat, eventInd * p_jk)
    aggregate(tildeN_jk, by = list(time = dat$time), FUN = sum)[, "x"]
}

k0 <- function(dat) {
    p_jk_xExp <- with(dat, p_jk * xExp)
    ## note that aggregate convert time into factor
    ## so the output will be sorted increasingly automatically
    uniDat <- aggregate(p_jk_xExp, by = list(time = dat$time), FUN = sum)
    ## vector of length # unique event time
    rev(cumsum(rev(uniDat$x)))
}

k1 <- function(parSeq, dat, xMat) {
    ## matrix of dimension # unique event time by # parameters
    sapply(parSeq, function(ind, dat, xMat) {
        p_jk_xExp_x <- with(dat, p_jk * xExp) * xMat[, ind]
        uniDat <- aggregate(p_jk_xExp_x, by = list(time = dat$time), FUN = sum)
        rev(cumsum(rev(uniDat$x)))
    }, dat = dat, xMat = xMat)
}

k2 <- function(parSeq, dat, xMat) {
    ind_grid <- expand.grid(parSeq, parSeq)
    ## matrix of dimension # unique event time by (# parameters) ^ 2
    mapply(function(ind1, ind2) {
        p_jk_xExp_x2 <- with(dat, p_jk * xExp) * xMat[, ind1] * xMat[, ind2]
        uniDat <- aggregate(p_jk_xExp_x2, by = list(time = dat$time), FUN = sum)
        rev(cumsum(rev(uniDat$x)))
    }, ind_grid[, 1L], ind_grid[, 2L])
}

dLbeta <- function(xMatDeltaN, k_0, k_1, delta_tildeN) {
    sum_jk <- colSums(xMatDeltaN)
    int_t <- colSums(k_1 / k_0 * delta_tildeN)
    sum_jk - int_t
}

d2Lbeta <- function(parSeq, k_0, k_1, k_2, delta_tildeN) {
    ## part 1
    nPar <- length(parSeq)
    mat1 <- k_2 / k_0 * delta_tildeN
    part1 <- matrix(colSums(mat1), nPar, nPar)

    ## part 2
    mat2 <- k_1 / k_0
    parGrid <- expand.grid(parSeq, parSeq)
    mat2plus <- mapply(function(ind1, ind2) {
        mat2[, ind1] * mat2[, ind2]
    }, parGrid[, 1L], parGrid[, 2L])
    part2 <- matrix(colSums(mat2plus * delta_tildeN), nPar, nPar)

    part2 - part1
}

coxEmStart <- function(beta, h0 = 1e-3, censorRate = 0.8, ..., nBeta_, dat_) {
    ## baseline hazard function: h0Vec
    ## subDat <- base::subset(dat, event == 0L)
    ## tempFit <- survival::survreg(survival::Surv(time, event) ~ X1 + X2,
    ##                              data = subDat, scale = 1)
    ## constant baseline hazard function
    ## lambda <- exp(- coef(tempFit)[1])

    ## initialize baseline hazard rate
    h0Vec <- ifelse(dat_$eventInd, h0, 0)

    ## initialize covariate coefficient: beta
    if (missing(beta)) {
        beta <- rep(0, nBeta_)
    } else if (length(beta) != nBeta_) {
        stop(paste("Number of starting values for coefficients of covariates",
                   "does not match with the specified formula."))
    } else {
        beta <- as.vector(beta, mode = "numeric")
    }

    ## mixture probability for each subject: piVec
    numTab <- table(dat_$ID)
    dat_$numRecord <- numTab[match(as.character(dat_$ID), names(numTab))]

    ## subject ID with censoring records
    cenID <- with(subset(dat_, ! eventInd), unique(ID))
    cenIdx <- as.integer(dat_$ID %in% cenID)

    ## for subjects with multiple records
    idx1 <- with(dat_, numRecord > 1L & eventInd)
    idx2 <- with(dat_, numRecord > 1L & (! eventInd))

    piVec <- rep(1L, NROW(dat_))
    piVec[idx1] <- (1 - censorRate * cenIdx[idx1]) /
        (dat_[idx1, "numRecord"] - cenIdx[idx1])
    piVec[idx2] <- censorRate

    ## return
    list(beta = beta, h0Vec = h0Vec, piVec = piVec)
}


coxEmControl <- function(gradtol = 1e-6, stepmax = 1e5,
                         steptol = 1e-6, iterlim = 1e2,
                         tolEm = 1e-6, iterlimEm = 1e3,
                         tolSem = 1e-3, iterlimSem = 1e2, ...) {
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
         tolSem = tolSem, iterlimSem = iterlimSem)
}
