################################################################################
### functions fitting Cox model with uncertain survival record data
### version controlled by git
################################################################################


## implementation of EM algorithm to Cox model
coxEm <- function(formula, data, subset, na.action, contrasts = NULL,
                  start = list(), control = list(), ...) {

    ## record the function call to return
    Call <- match.call()

    ## arguments check
    if (missing(formula)) {
        stop("Argument 'formula' is required.")
    }
    if (missing(data)) {
        data <- environment(formula)
    }
    if (! with(data, inherits(eval(Call[[2]][[2]]), "Surve"))) {
        stop("Response in formula must be a 'surve' object.")
    }

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
    if ((nBeta <- ncol(mm) - 1L) <= 0) {
        stop("Covariates must be specified in formula.")
    }

    ## covariates' names
    covar_names <- colnames(mm)[-1]

    ## data
    dat <- as.data.frame(cbind(mf[, 1][, 1:3], mm[, -1]))
    colnames(dat) <- c("ID", "time", "event", covar_names)
    outData <- dat
    nObs <- nrow(dat)
    nBeta <- ncol(dat) - 3L
    xMat <- as.matrix(dat[, - (1:3)])

    ## define some variables for ease of computing
    dat$eventInd <- dat$event == 1L
    dat$orderInc <- order(dat$time, decreasing = FALSE)
    dat$orderDec <- order(dat$time, decreasing = TRUE)

    ## 'control' for 'nlm'
    control <- do.call("coxEmControl", control)

    ## start' values for 'nlm'
    startList <- c(start, list(nBeta = nBeta, dat = dat))
    start <- do.call("coxEmStart", startList)
    betaHat <- start$beta
    h0Vec <- start$h0Vec

    ## pkjVec <- NULL
    piVec <- start$piVec
    iter <- tol <- lMin <- 1L
    logL <- rep(NA, control$iterlimEm)
    while (iter < control$iterlimEm && tol > control$tolEm) {
        oneFit <- do1mStep(betaHat = betaHat, h0Vec = h0Vec, piVec = piVec,
                           dat = dat, xMat = xMat, control = control)
        piVec <- pkjVec <- oneFit$pkjVec
        h0Vec <- oneFit$h0Vec
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
    est_beta[, 4] <- est_beta[, 1]/est_beta[, 3]
    est_beta[, 5] <- 2 * stats::pnorm(- abs(est_beta[, 4]))

    ## output: na.action
    if (is.null(attr(mf, "na.action"))) {
        na.action <- options("na.action")[[1]]
    } else {
        na.action <- paste("na", class(attr(mf, "na.action")), sep = ".")
    }
    ## output: contrasts
    if (is.null(contrasts)) {
        contrasts <- list(contrasts = NULL)
    } else {
        contrasts <- attr(mm, "contrasts")
    }
    ## results to return
    results <- methods::new("coxEm",
                            call = Call,
                            formula = formula,
                            data = outData,
                            nObs = nObs,
                            estimates = list(beta = est_beta,
                                             pkjVec = pkjVec,
                                             h0Vec = h0Vec),
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
do1mStep <- function(betaHat, h0Vec, piVec, dat, xMat, control) {

    dat$betaX <- xMat %*% betaHat
    dat$xExp <- pmax(exp(dat$betaX), .Machine$double.eps)
    h0Vec <- H0Vec <- ifelse(dat$eventInd, h0Vec, 0)
    H0Vec[dat$orderInc] <- cumsum(h0Vec[dat$orderInc])
    hVec <- h0Vec * dat$xExp
    HVec <- H0Vec * dat$xExp
    sVec <- pmax(exp(- HVec), .Machine$double.eps)

    ## compute p_kj for each subject
    p_kj_numer <- ifelse(dat$eventInd,
                         pmax(piVec * hVec * sVec, .Machine$double.eps), 1)
    tempDat <- cbind(dat, hVec, HVec, sVec, p_kj_numer)
    p_kj_denom_dat <- aggregate(p_kj_numer ~ ID, data = tempDat, FUN = sum)
    colnames(p_kj_denom_dat)[2] <- "p_kj_denom"
    outDat <- merge(tempDat, p_kj_denom_dat, by = "ID")
    outDat$p_kj <- with(outDat, p_kj_numer / p_kj_denom)

    ## likelihood function component for each subject
    outDat$subL <- with(outDat, ifelse(eventInd, p_kj_denom, sVec))
    ## likelihood function
    logL <- with(outDat, sum(log(subL[! duplicated(ID)])))

    ## update beta
    betaEst <- stats::nlm(logLbeta, p = betaHat, dat = outDat, xMat = xMat,
                          hessian = TRUE, check.analyticals = TRUE,
                          gradtol = control$gradtol, stepmax = control$stepmax,
                          steptol = control$steptol, iterlim = control$iterlim)

    ## compute with updated beta estimates
    outDat$betaX <- xMat %*% betaEst$estimate
    outDat$xExp <- pmax(exp(outDat$betaX), .Machine$double.eps)
    ## update h0_kj with previous (or initial) estimates of beta
    outDat$h0_kj <- h0kj(outDat)

    ## return
    list(betaEst = betaEst, h0Vec = outDat$h0_kj,
         pkjVec = outDat$p_kj, logL = logL)
}

h0kj <- function(dat) {
    sortDat <- dat[dat$orderDec, ]
    dat$h0_kj_denom <- NA
    dat$h0_kj_denom[dat$orderDec] <- with(sortDat, cumsum(p_kj * xExp))
    with(dat, p_kj * event / h0_kj_denom)
}

dLbeta <- function(par, dat, xMat) {
    dLbeta1 <- function(ind, dat, xMat) {
        xTime <- dat$time
        sortDat <- dat[dat$orderDec, ]
        denom <- with(sortDat, cumsum(p_kj * xExp))
        numer <- cumsum(sortDat$p_kj * sortDat$xExp * xMat[dat$orderDec, ind])
        tempPart <- rep(NA, length(xTime))
        tempPart[dat$orderDec] <- numer / denom
        res <- (dat$p_kj * (xMat[, ind] - tempPart))[dat$eventInd]
        sum(res)
    }
    sapply(seq_along(par), dLbeta1, dat = dat, xMat = xMat)
}

d2Lbeta <- function(par, dat, xMat) {
    d2Lbeta1 <- function(ind1, ind2, dat, xMat) {
        xTime <- dat$time
        sortDat <- dat[dat$orderDec, ]
        xVec1 <- xMat[dat$orderDec, ind1]
        xVec2 <- xMat[dat$orderDec, ind2]
        denom <- with(sortDat, cumsum(p_kj * xExp))
        numer1 <- with(sortDat, cumsum(p_kj * xExp * xVec1)) *
            with(sortDat, cumsum(p_kj * xExp * xVec2))
        numer2 <- with(sortDat, cumsum(p_kj * xExp * xVec1 * xVec2))
        tempPart <- rep(NA, length(xTime))
        tempPart[dat$orderDec] <- numer1 / denom ^ 2 - numer2 / denom
        res <- (dat$p_kj * tempPart)[dat$eventInd]
        sum(res)
    }
    hesMat <- matrix(NA, nrow = length(par), ncol = length(par))
    for (d1 in seq_along(par)) {
        for (d2 in seq_along(par)) {
            hesMat[d1, d2] <- d2Lbeta1(d1, d2, dat = dat, xMat = xMat)
        }
    }
    hesMat
}

logLbeta <- function(par, dat, xMat) {

    dat$betaX <- xMat %*% par
    dat$xExp <- pmax(exp(dat$betaX), .Machine$double.eps)
    h0_kj <- h0kj(dat)
    tempPart <- ifelse(dat$eventInd & (h0_kj != 0L),
                       dat$p_kj * (log(h0_kj) + dat$betaX), 0L)

    ## penalty term to avoid solution with xExp being Inf
    penal <- any(dat$xExp == Inf) * 1e20

    negLogL <- - sum(tempPart) + penal
    ## gradient
    gradLogL <- dLbeta(par = par, dat = dat, xMat = xMat)
    attr(negLogL, "gradient") <- - gradLogL
    ## hessian
    hesMat <- d2Lbeta(par = par, dat = dat, xMat = xMat)
    attr(negLogL, "hessian") <- hesMat * (- 1)
    ## return
    negLogL
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
    if (!is.numeric(gradtol) || gradtol <= 0) {
        stop("value of 'gradtol' must be > 0")
    }
    if (!is.numeric(stepmax) || stepmax <= 0) {
        stop("value of 'stepmax' must be > 0")
    }
    if (!is.numeric(steptol) || steptol <= 0) {
        stop("value of 'steptol' must be > 0")
    }
    if (!is.numeric(iterlim) || iterlim <= 0) {
        stop("maximum number of iterations must be > 0")
    }
    if (!is.numeric(tolEm) || tolEm <= 0) {
        stop("value of 'tolEm' must be > 0")
    }
    if (!is.numeric(iterlimEm) || iterlimEm <= 0) {
        stop("maximum number of iterations for EM must be > 0")
    }
    ## return
    list(gradtol = gradtol, stepmax = stepmax,
         steptol = steptol, iterlim = iterlim,
         tolEm = tolEm, iterlimEm = iterlimEm)
}




