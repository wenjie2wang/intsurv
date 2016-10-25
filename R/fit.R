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
##'                    censorTime = 3.48, lambda = 1e-2, rho = 2, mixture = 0.5)
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
    startList <- c(start, list(nBeta = nBeta, dat = incDat))
    start <- do.call("coxEmStart", startList)
    betaHat <- start$beta

    ## take care of possible ties
    h0Dat <- aggregate(start$h0Vec, list(time = dat$time), FUN = sum)
    colnames(h0Dat)[2L] <- "h0Vec"

    piVec <- start$piVec[orderInc]
    iter <- tol <- lMin <- 1L
    logL <- rep(NA, control$iterlimEm)
    while (iter < control$iterlimEm && tol > control$tolEm) {
        oneFit <- oneEMstep(betaHat = betaHat, h0Dat = h0Dat, piVec = piVec,
                            dat = incDat, xMat = xMat, control = control)
        piVec <- pjkVec <- oneFit$pjkVec
        h0Dat$h0Vec <- oneFit$h0Vec
        betaEst <- oneFit$betaEst
        betaHatNew <- betaEst$estimate
        ## log likehood
        logL[iter] <- oneFit$logL
        iter <- iter + 1
        tol <- sqrt(sum((betaHatNew - betaHat) ^ 2)) / sqrt(sum(betaHat^2))
        ## update beta estimates
        betaHat <- betaHatNew
    }

    ## log likelihood
    logL <- stats::na.omit(logL)
    attr(logL, "na.action") <- attr(logL, "class") <- NULL

    ## estimates for beta
    est_beta <- as.matrix(betaHat)
    est_beta <- matrix(NA, nrow = nBeta, ncol = 5)
    colnames(est_beta) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
    rownames(est_beta) <- covar_names
    se_vec <- sqrt(diag(solve(betaEst$hessian)))
    est_beta[, 1] <- betaEst$estimate[1 : nBeta]
    est_beta[, 2] <- exp(est_beta[, 1])
    est_beta[, 3] <- se_vec[1:nBeta]
    est_beta[, 4] <- est_beta[, 1] / est_beta[, 3]
    est_beta[, 5] <- 2 * stats::pnorm(- abs(est_beta[, 4]))

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
                                             pjkVec = pjkVec,
                                             h0Dat = h0Dat),
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
## perform one step of EM algorithm
oneEMstep <- function(betaHat, h0Dat, piVec, dat, xMat, control) {

    dat$betaX <- as.vector(xMat %*% betaHat)
    dat$xExp <- pmax(exp(dat$betaX), .Machine$double.eps)
    ## the order is preserved
    tmpDat <- merge(dat, h0Dat, by = "time", sort = FALSE)

    ## update baseline hazard rate
    tmpDat$h0Vec <- with(tmpDat, ifelse(eventInd, h0Vec, 0))
    tmpDat$H0Vec <- cumsum(tmpDat$h0Vec)
    tmpDat$hVec <- with(tmpDat, h0Vec * xExp)
    tmpDat$HVec <- with(tmpDat, H0Vec * xExp)
    tmpDat$sVec <- pmax(exp(- tmpDat$HVec), .Machine$double.eps)

    ## compute p_jk for each subject
    tmpDat$p_jk_numer <- with(tmpDat, ifelse(eventInd, piVec * hVec * sVec,
                                             piVec * sVec))
    tmpDat$p_jk_numer <- pmax(tmpDat$p_jk_numer, .Machine$double.eps)
    p_jk_denom_dat <- aggregate(p_jk_numer ~ ID, data = tmpDat, FUN = sum)
    colnames(p_jk_denom_dat)[2] <- "p_jk_denom"

    ## record originally sorted order
    tmpDat$idx <- seq_len(nrow(tmpDat))
    outDat <- merge(tmpDat, p_jk_denom_dat, by = "ID", sort = FALSE)
    outDat$p_jk <- with(outDat, p_jk_numer / p_jk_denom)
    ## recover the order for piVec estimates
    outDat <- outDat[order(outDat$idx), ]

    ## update design matrix
    parSeq <- seq_along(betaHat)
    xMat <- as.matrix(outDat[, 3L + parSeq])

    ## log-likelihood function under observed data
    logL <- with(outDat, sum(log(p_jk_denom[! duplicated(ID)])))

    ## update beta
    betaEst <- stats::nlm(logLbeta, p = betaHat, dat = outDat, xMat = xMat,
                          hessian = TRUE, check.analyticals = TRUE,
                          gradtol = control$gradtol, stepmax = control$stepmax,
                          steptol = control$steptol, iterlim = control$iterlim)

    ## compute with updated beta estimates
    outDat$betaX <- as.vector(xMat %*% betaEst$estimate)
    outDat$xExp <- pmax(exp(outDat$betaX), .Machine$double.eps)
    ## update h0_jk with previous (or initial) estimates of beta
    h0Vec <- h0t(outDat)

    ## return
    list(betaEst = betaEst, h0Vec = h0Vec, pjkVec = outDat$p_jk, logL = logL)
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

coxEmStart <- function(beta, ..., nBeta, dat) {
    ## baseline hazard function: h0Vec
    ## subDat <- base::subset(dat, event == 0L)
    ## tempFit <- survival::survreg(survival::Surv(time, event) ~ X1 + X2,
    ##                              data = subDat, scale = 1)
    ## constant baseline hazard function
    ## lambda <- exp(- coef(tempFit)[1])
    lambda <- 1e-3
    h0Vec <- rep(lambda, nrow(dat))

    ## covariate coefficient: beta
    if (missing(beta)) {
        beta <- rep(0, nBeta)
    } else if (length(beta) != nBeta) {
        stop(paste("Number of starting values for coefficients of covariates",
                   "does not match with the specified formula."))
    } else {
        beta <- as.vector(beta, mode = "numeric")
    }

    ## mixture probability for each subject: piVec
    numTab <- table(dat$ID)
    piVec <- rep(1 / numTab, numTab)

    ## return
    list(beta = beta, h0Vec = h0Vec, piVec = piVec)
}


coxEmControl <- function(gradtol = 1e-6, stepmax = 1e5,
                         steptol = 1e-6, iterlim = 1e2,
                         tolEm = 1e-6, iterlimEm = 1e3, ...) {
    ## controls for function stats::nlm
    if (!is.numeric(gradtol) || gradtol <= 0)
        stop("value of 'gradtol' must be > 0")
    if (!is.numeric(stepmax) || stepmax <= 0)
        stop("value of 'stepmax' must be > 0")
    if (!is.numeric(steptol) || steptol <= 0)
        stop("value of 'steptol' must be > 0")
    if (!is.numeric(iterlim) || iterlim <= 0)
        stop("maximum number of iterations must be > 0")
    if (!is.numeric(tolEm) || tolEm <= 0)
        stop("value of 'tolEm' must be > 0")
    if (!is.numeric(iterlimEm) || iterlimEm <= 0)
        stop("maximum number of iterations for EM must be > 0")
    ## return
    list(gradtol = gradtol, stepmax = stepmax,
         steptol = steptol, iterlim = iterlim,
         tolEm = tolEm, iterlimEm = iterlimEm)
}
