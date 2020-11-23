##
## intsurv: Integrative Survival Models
## Copyright (C) 2017-2020  Wenjie Wang <wang@wwenjie.org>
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


##' Integrative Cox Model for Uncertain Event Times
##'
##' Fit an integrative Cox model proposed by Wang et al. (2020) for
##' right-censored survival data with uncertain event times due to imperfect
##' data integration.
##'
##' @param formula \code{Survi} object specifying the covariates and response
##'     variable in the model, such as \code{Survi(ID, time, event) ~ x1 + x2}.
##' @param data An optional data frame, list, or environment that contains the
##'     covariates and response variables included in the model. If not found in
##'     data, the variables are taken from \code{environment(formula)}, usually
##'     the environment from which this function is called.
##' @param subset An optional logical vector specifying a subset of observations
##'     to be used in the fitting process.
##' @param na.action An optional function that indicates what should the
##'     procedure do if the data contains \code{NA}s.  The default is set by the
##'     na.action setting of \code{\link[base]{options}}.  The "factory-fresh"
##'     default is \code{\link[stats]{na.omit}}.  Other possible values include
##'     \code{\link[stats]{na.fail}}, \code{\link[stats]{na.exclude}}, and
##'     \code{\link[stats]{na.pass}}.  \code{help(na.fail)} for details.
##' @param contrasts An optional list, whose entries are values (numeric
##'     matrices or character strings naming functions) to be used as
##'     replacement values for the contrasts replacement function and whose
##'     names are the names of columns of data containing factors.  See
##'     \code{contrasts.arg} of \code{\link[stats]{model.matrix.default}} for
##'     details.
##' @param start A list returned by function \code{\link{iCoxph.start}}
##'     specifying starting values of the parameters to be estimated in the
##'     model. Please refer to the arguments of \code{\link{iCoxph.start}} for
##'     the available parameters.
##' @param control A list returned by function \code{\link{iCoxph.control}}
##'     specifying control parameters for the model estimation procedure.
##'     Please refer to the arguments of \code{\link{iCoxph.control}} for the
##'     available parameters.
##' @param ... Other arguments for future usage.  A warning will be thrown if
##'     any invalid argument is specified.
##'
##' @return
##' An \code{\link{iCoxph-class}} object, whose slots include
##' \itemize{
##'     \item \code{call}: Function call.
##'     \item \code{formula}: Formula used in the model fitting.
##'     \item \code{data}: A processed data frame used for model fitting.
##'     \item \code{nObs}: Number of observation.
##'     \item \code{estimates}:
##'         \itemize{
##'             \item \code{beta}: Coefficient estimates.
##'             \item \code{pi}: Estimated parameters in prior multinomial
##'                 distribution indicating the probabilities of corresponding
##'                 record being true.
##'             \item \code{baseline}: A data frame that contains estimated
##'                 baseline hazard function with columns named \code{time} and
##'                 \code{h0}.
##'         }
##'     \item \code{start}: The initial guesses \code{beta_mat} and
##'         \code{pi_mat} specified for the parameters to be estimated,
##'         and the set of initial guess \code{beta_start} and \code{pi_start}
##'         that resulted in the largest objective function, i.e.,
##'         the observed-data likelihood function.
##'     \item \code{control}: The control list specified for model fitting.
##'     \item \code{na.action}: The procedure specified to deal with
##'         missing values in the covariate.
##'     \item \code{xlevels}: A list that records the levels in
##'         each factor variable.
##'     \item \code{contrasts}: Contrasts specified and used for each
##'         factor variable.
##'     \item \code{convergeCode}: \code{code} returned by function
##'         \code{\link[stats]{nlm}}, which is an integer indicating why the
##'         optimization process terminated. \code{help(nlm)} for details.
##'     \item \code{logL}: A numeric vector containing the observed-data
##'         log-likelihood over iterations.
##' }
##'
##' @references
##'
##' Wang, W., Aseltine, R. H., Chen, K., & Yan, J. (2020). Integrative Survival
##' Analysis with Uncertain Event Times in Application to A Suicide Risk
##' Study. \emph{Annals of Applied Statistics}, 14(1), 51--73.
##'
##' @example inst/examples/iCoxph.R
##'
##' @seealso
##' \code{\link{iCoxph.start}} and \code{\link{iCoxph.control}}, respectively,
##' for starting and controlling iCoxph fitting;
##' \code{\link{summary,iCoxph-method}} for summary of fitted model;
##' \code{\link{coef,iCoxph-method}} for estimated covariate coefficients;
##' \code{\link{bootSe}} for SE estimates from bootstrap methods.
##'
##' @importFrom stats na.fail na.omit na.exclude na.pass .getXlevels
##'     model.extract model.frame model.matrix nlm pnorm
##'
##' @export
iCoxph <- function(formula, data, subset, na.action, contrasts = NULL,
                   start = iCoxph.start(), control = iCoxph.control(), ...)
{
    ## record the function call to return
    Call <- match.call()

    ## warning on `...`
    warn_dots(...)

    ## arguments check
    if (missing(formula))
        stop("Argument 'formula' is required.")
    if (missing(data))
        data <- environment(formula)
    dat0 <- with(data, eval(formula[[2L]]))
    if (! is_Survi(dat0))
        stop("The formula response must be a 'Survi' object.")
    Call$formula <- formula

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

    ## prepare start' values
    if (! is_iCoxph.start(start)) {
        stop("The 'start' must be specified through ",
             "function 'iCoxph.start()'")
    }
    startList <- c(start, list(nBeta_ = nBeta, dat_ = dat))
    start <- unclass(do.call("iCoxph_start", startList))

    ## prepare 'control' parameters
    if (! is_iCoxph.control(control)) {
        stop("The 'control' must be specified through ",
             "function 'iCoxph.control()'")
    }
    control <- unclass(control)

    ## sort by time and ID
    incDat <- dat[(orderInc <- with(dat, order(time, ID))), ]
    ## indicator for subjects having multiple record
    incDat$dupIdx <- with(incDat, ID %in% unique(ID[duplicated(ID)]))
    ## define some variables for ease of computing
    incDat$eventIdx <- incDat$event == 1L
    dupVec <- duplicated(incDat$time)
    tied <- any(dupVec)
    incDat$firstIdx <- ! dupVec
    h0Dat <- h_cDat <- data.frame(time = incDat$time[incDat$firstIdx])
    xMat <- as.matrix(incDat[, 4L : (3L + nBeta)])
    piMat <- start$pi_mat[orderInc, , drop = FALSE]

    ## initialize log-likelihood
    logL_max0 <- - Inf
    n_beta_start <- ncol(start$beta_mat)
    n_pi_start <- ncol(piMat)
    logL_max_vec <- rep(NA, n_beta_start * n_pi_start)
    beta_est_mat <- matrix(NA, nrow = n_beta_start * n_pi_start, ncol = nBeta)

    for (oneBeta in seq_len(n_beta_start)) {
        for (onePi in seq_len(n_pi_start)) {
            ## index of starts
            idx_start <- 1
            incDat$piVec <- piVec <- piMat[, onePi]
            ## trace the log-likelihood for observed data
            logL <- rep(NA, control$ecm_max_iter)
            ## trace beta estimates from each iteration of ECM
            betaMat <- matrix(NA, nrow = control$ecm_max_iter + 1L,
                               ncol = nBeta)
            betaMat[1L, ] <- start$beta_mat[, oneBeta]
            tol_update <- sqrt(control$tol_beta)

            for (iter in seq_len(control$ecm_max_iter)) {
                oneFit <- ic_oneECMstep(betaHat = betaMat[iter, ],
                                        h0Dat = h0Dat,
                                        h_cDat = h_cDat,
                                        dat = incDat,
                                        xMat = xMat,
                                        tied = tied,
                                        control = control)
                ## log likehood
                logL[iter] <- oneFit$logL

                ## update pi estimates
                incDat$piVec <- oneFit$piVec

                ## update beta estimates
                betaEst <- oneFit$betaEst
                betaMat[iter + 1L, ] <- betaEst$estimate
                tol_beta <- L2norm(betaMat[iter + 1L, ] - betaMat[iter, ]) /
                    L2norm(betaMat[iter + 1L, ] + betaMat[iter, ])
                tol_pi <- L2norm(oneFit$piVec - incDat$piVec) /
                    L2norm(oneFit$piVec + incDat$piVec)

                betaHat <- betaEst$estimate
                if (tol_beta < control$tol_beta &&
                    tol_pi < control$tol_pi) {
                    break
                } else if (iter == control$ecm_max_iter) {
                    warning("Reached the maximum number of ECM iterations.")
                }
            }

            ## keep the one fit maximizing observed log likelihood
            logL <- rmNA(logL)
            iter <- (oneBeta - 1) * n_pi_start + onePi
            n_step <- length(logL)
            logL_max <- logL_max_vec[idx_start] <- logL[n_step]
            beta_est_mat[iter, ] <- betaHat

            if (logL_max > logL_max0) {
                logL_max0 <- logL_max
                logL0 <- logL
                oneFit0 <- oneFit
                betaMat0 <- betaMat
                piVec0 <- piVec
                start$beta_start <- betaMat[1L, ]
            }

            ## update index of start
            idx_start <- idx_start + 1
        }
    }

    ## clean-up NA's
    betaMat0 <- rmNA(betaMat0)
    ## prepare for outputs
    piEst <- oneFit0$piVec[(reOrderIdx <- order(orderInc))]
    start$pi_start <- piVec0[reOrderIdx]
    ## start$logL_end <- logL_max_vec
    ## start$beta_mat_end <- t(beta_est_mat)
    ## update results
    h0Dat$h0Vec <- oneFit0$h0Vec
    colnames(h0Dat) <- c("time", "h0")
    h_cDat$h_cVec <- oneFit0$h_cVec
    incDat$xExp <- oneFit0$xExp
    incDat$piVec <- oneFit0$piVec
    betaEst <- oneFit0$betaEst
    betaHat <- as.numeric(betaEst$estimate)

    ## estimates for beta
    est_beta <- matrix(NA, nrow = nBeta, ncol = 5L)
    colnames(est_beta) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
    rownames(est_beta) <- covar_names
    est_beta[, "coef"] <- betaHat
    est_beta[, "exp(coef)"] <- exp(est_beta[, "coef"])
    est_beta[, "se(coef)"] <- NA_real_
    est_beta[, "z"] <- est_beta[, "coef"] / est_beta[, "se(coef)"]
    est_beta[, "Pr(>|z|)"] <- 2 * stats::pnorm(- abs(est_beta[, "z"]))

    ## output: processed data frame
    dat$ID <- factor(dat$ID, levels = unique(dat$ID),
                     labels = unique(dat0@ID))
    dat$ID <- as.character(dat$ID)
    varNames0 <- sapply(seq_len(3L), function(a) {
        as.character(formula[[2L]][[a + 1]])
    })
    colnames(dat)[seq_len(3L)] <- varNames0

    ## output: na.action
    na.action <- if (is.null(attr(mf, "na.action")))
                     options("na.action")[[1L]]
                 else
                     paste0("na.", class(attr(mf, "na.action")))

    ## results to return
    methods::new("iCoxph",
                 call = Call,
                 formula = formula,
                 nObs = nObs,
                 data = dat,
                 estimates = list(beta = est_beta,
                                  pi = piEst,
                                  baseline = h0Dat),
                 start = start,
                 control = control,
                 na.action = na.action,
                 xlevels = .getXlevels(mt, mf),
                 contrasts = list(attr(mm, "contrasts")),
                 convergeCode = betaEst$code,
                 logL = logL0)
}


##' Auxiliary for Controlling iCoxph Fitting
##'
##' Auxiliary function for \code{\link{iCoxph}} that enable users to specify the
##' control parameters of the model estimation procedure.  Internally, the
##' arguments \code{cm_gradtol}, \code{cm_stepmax}, \code{cm_steptol}, and
##' \code{cm_max_iter} are passed to function \code{\link[stats]{nlm}} as
##' \code{gradtol}, \code{stepmax}, \code{steptol}, and \code{iterlim},
##' respectively.
##'
##' @param tol_beta A positive value specifying the tolerance that concludes the
##'     convergence of the covariate coefficient estimates. The tolerance is
##'     compared with the relative change between the estimates from two
##'     consecutive iterations that is measured by ratio of the L2-norm of their
##'     difference to the sum of their L2-norm.  The default value is
##'     \code{1e-6}.
##' @param tol_pi A positive value specifying the tolerance that concludes the
##'     convergence of the probability estimates of uncertain records being
##'     true.  The tolerance is compared with the relative change between the
##'     estimates from two consecutive iterations measured by ratio of L2-norm
##'     of their difference to the L2-norm of their sum.  The default value is
##'     \code{1e-8}.
##' @param cm_gradtol A positive scalar giving the tolerance at which the scaled
##'     gradient is considered close enough to zero to terminate CM steps.  The
##'     default value is \code{1e-6}.
##' @param cm_stepmax A positive scalar that gives the maximum allowable scaled
##'     step length in CM steps.  The default value is \code{1e2}.
##' @param cm_steptol A positive scalar providing the minimum allowable relative
##'     step length in CM step.  The default value is \code{1e-6}.
##' @param cm_max_iter A positive integer specifying the maximum number of
##'     iterations to be performed before each CM step is terminated. The
##'     default value is \code{1e2}.
##' @param ecm_max_iter A positive integer specifying the maximum number of
##'     iterations to be performed before the ECM algorithm is terminated. The
##'     default value is \code{2e2}.
##' @param ... Other arguments for future usage.  A warning will be thrown if
##'     any invalid argument is specified.
##'
##' @return A list of class \code{intsurv-iCoxph.control} containing all
##'     specified control parameters.
##'
##' @examples
##' ## See examples of function 'iCoxph'.
##'
##' @seealso
##' \code{\link{iCoxph}} for fitting integrative Cox model.
##'
##' @export
iCoxph.control <- function(tol_beta = 1e-6,
                           tol_pi = 1e-8,
                           cm_gradtol = 1e-6,
                           cm_stepmax = 1e2,
                           cm_steptol = 1e-6,
                           cm_max_iter = 1e2,
                           ecm_max_iter = 2e2,
                           ...)
{
    ## controls for function stats::nlm
    if (! is.numeric(cm_gradtol) || cm_gradtol <= 0)
        stop("The value of 'cm_gradtol' must be > 0.", call. = FALSE)
    if (! is.numeric(cm_stepmax) || cm_stepmax <= 0)
        stop("The value of 'cm_stepmax' must be > 0.", call. = FALSE)
    if (! is.numeric(cm_steptol) || cm_steptol <= 0)
        stop("The value of 'cm_steptol' must be > 0.", call. = FALSE)
    if (! is.numeric(cm_max_iter) || cm_max_iter <= 0)
        stop("The maximum number of iterations in the CM step must be > 0.",
             call. = FALSE)
    ## determining convergence of ECM
    if (! is.numeric(tol_beta) || tol_beta <= 0)
        stop("The value of 'tol_beta' must be positive.", call. = FALSE)
    if (! is.numeric(tol_pi) || tol_pi <= 0)
        stop("The value of 'tol_pi' must be positive.", call. = FALSE)
    if (! is.numeric(ecm_max_iter) || ecm_max_iter <= 0)
        stop("The maximum number of ECM iterations must be positive.",
             call. = FALSE)
    ## throw warning if anything is captured by `...`
    warn_dots(...)
    ## parameters in a list
    out <- list(tol_beta = tol_beta,
                tol_pi = tol_pi,
                cm_gradtol = cm_gradtol,
                cm_stepmax = cm_stepmax,
                cm_steptol = cm_steptol,
                cm_max_iter = cm_max_iter,
                ecm_max_iter = ecm_max_iter)
    ## add class for ease of checking
    class(out) <- "iCoxph.control"
    ## return
    out
}


##' Auxiliary for Starting iCoxph Fitting
##'
##' Auxiliary function for \code{\link{iCoxph}} that enable users
##' to specify the starting values of the model estimation procedure.
##'
##' @param beta_vec A numeric vector for starting values of coefficient
##'     estimates. The default values are the coefficient estimates from the
##'     regular Cox model only fitting on records without uncertainty.  If
##'     censoring rate among subjects having unique certain records is extremely
##'     high (> 99%) or such as perfect seperation problem between event
##'     indicator and one predictor, the starting values will be reset to be all
##'     zeros.
##' @param beta_mat A numeric matrix that consists of additional starting values
##'     of coefficient estimates in columns.  The default value is \code{NULL}.
##' @param methods A character vector specifying the initialization methods for
##'     probabilities of uncertain records being true.  The available methods
##'     are \code{"nearest_hazard"} for initializing baseline hazard by nearest
##'     (left) neighbor, and \code{"unit_hazard"} for initializing unit baseline
##'     hazard.  Partial matching on method names is supported for ease of
##'     typing.  By default, both methods are used.  See Wang et al. (2020) for
##'     details of the initialization methods.
##' @param ... Other arguments for future usage.  A warning will be thrown if
##'     any invalid argument is specified.
##'
##' @return A list of class \code{intsurv-iCoxph.start} containing all specified
##'     starting values of the parameters to be estimated from the model.
##'
##' @examples
##' ## See examples of function 'iCoxph'.
##'
##' @seealso
##' \code{\link{iCoxph}} for fitting integrative Cox model.
##'
##' @export
iCoxph.start <- function(beta_vec = NULL,
                         beta_mat = NULL,
                         methods = c("nearest_hazard", "unit_hazard"),
                         ...)
{
    ## match methods
    methods <- match.arg(methods, several.ok = TRUE)
    ## throw out warning if any invalid argument is captured by `...`
    warn_dots(...)
    ## parameters in a list
    out <- list(beta_vec = beta_vec,
                beta_mat = beta_mat,
                methods = methods)
    ## add class for ease of checking
    class(out) <- "iCoxph.start"
    ## return
    out
}


### internal functions =========================================================
## perform one step of ECM algorithm
ic_oneECMstep <- function(betaHat, h0Dat, h_cDat, dat, xMat, tied, control)
{
    ## update results involving beta estimates
    dat$betaX <- as.numeric(xMat %*% betaHat)
    dat$xExp <- exp(dat$betaX)
    dat$xExp <- ifelse(is.infinite(dat$xExp), 1e50, dat$xExp)

    ## use initial value of pi
    dat$p_jk <- dat$piVec
    h0Dat$h0Vec <- ic_h0t(dat, tied = tied)
    h_cDat$h_cVec <- ic_h_c(dat, tied = tied)
    ## help converge more quickly
    dat$p_jk <- ifelse(dat$p_jk < 1e-3, 0, ifelse(dat$p_jk > 1 - 1e-3,
                                                  1, dat$p_jk))

    ## update baseline hazard rate of event times
    h0Dat$H0Vec <- cumsum(h0Dat$h0Vec)
    ## update baseline hazard rate of censoring times
    h_cDat$H_cVec <- cumsum(h_cDat$h_cVec)

    if (tied) {
        time_idx <- match(dat$time, h0Dat$time)
        dat$h0Vec <- with(dat, ifelse(eventIdx, h0Dat$h0Vec[time_idx], 0))
        dat$H0Vec <- h0Dat$H0Vec[time_idx]
        dat$h_cVec <- with(dat, ifelse(eventIdx, 0, h_cDat$h_cVec[time_idx]))
        dat$H_cVec <- h_cDat$H_cVec[time_idx]
    } else {
        dat$h0Vec <- with(dat, ifelse(eventIdx, h0Dat$h0Vec, 0))
        dat$H0Vec <- h0Dat$H0Vec
        dat$h_cVec <- with(dat, ifelse(eventIdx, 0, h_cDat$h_cVec))
        dat$H_cVec <- h_cDat$H_cVec
    }

    dat$hVec <- with(dat, h0Vec * xExp)
    dat$HVec <- with(dat, H0Vec * xExp)
    dat$sVec <- exp(- dat$HVec)
    dat$G_cVec <- exp(- dat$H_cVec)

    ## compute p_jk for each subject
    ## for observed log-likelihood function
    dat$w_jk <- with(dat, ifelse(eventIdx,
                                 piVec * hVec * sVec * G_cVec,
                                 piVec * sVec * h_cVec * G_cVec))

    dat$p_jk_denom <- with(dat, aggregateSum(w_jk, ID, simplify = FALSE))
    dat$p_jk <- with(dat, ifelse(dupIdx, w_jk / p_jk_denom, 1))

    ## update h_0(t)
    h0Dat$h0Vec <- ic_h0t(dat, tied = tied)

    ## update h_c(t)
    h_cDat$h_cVec <- ic_h_c(dat, tied = tied)

    ## update beta
    betaEst <- stats::nlm(ic_logLbeta, p = betaHat,
                          dat = dat, xMat = xMat,
                          tied = tied,
                          hessian = FALSE,
                          check.analyticals = FALSE,
                          gradtol = control$cm_gradtol,
                          stepmax = control$cm_stepmax,
                          steptol = control$cm_steptol,
                          iterlim = control$cm_max_iter)

    ## log-likelihood function under observed data
    logL <- sum(log(dat$p_jk_denom))

    ## update h0_jk and h_c_jk with previous (or initial) estimates of beta
    ## h0Vec <- ic_h0t(dat, tied)
    ## h_cVec <- ic_h_c(dat, tied)
    list(betaEst = betaEst,
         h0Vec = h0Dat$h0Vec,
         h_cVec = h_cDat$h_cVec,
         piVec = dat$p_jk,
         logL = logL,
         xExp = dat$xExp)
}


## profile log-likelihood function of beta
ic_logLbeta <- function(param, dat, xMat, tied)
{
    ## update retults depends on beta estimates, param
    dat$betaX <- as.vector(xMat %*% param)
    dat$xExp <- exp(dat$betaX)
    dat$xExp <- ifelse(is.infinite(dat$xExp), 1e50, dat$xExp)

    ## prepare intermediate results for later computation
    parSeq <- seq_along(param)
    xMatDeltaN <- xMat[dat$eventIdx, ] * dat[dat$eventIdx, "p_jk"]
    delta_tildeN <- ic_deltaTildeN(dat, tied)
    betaXdeltaN <- with(dat, eventIdx * betaX * p_jk)

    k_0 <- ic_k0(dat, tied)
    k_1 <- ic_k1(parSeq, dat, xMat, tied)
    k_2 <- ic_k2(parSeq, dat, xMat, tied)

    ## profile log-likelihood of beta in EM
    pell <- sum(betaXdeltaN) - sum(ifelse(delta_tildeN > 0,
                                          delta_tildeN * log(k_0), 0))

    ## penalty term to avoid solution with xExp being Inf
    penal_inf <- any(is.infinite(dat$xExp)) * 1e20
    negLogL <- - pell + penal_inf

    ## gradient
    gradLogL <- ic_dLbeta(xMatDeltaN, k_0, k_1, delta_tildeN)
    attr(negLogL, "gradient") <- - gradLogL
    ## hessian
    hesMat <- ic_d2Lbeta(parSeq, k_0, k_1, k_2, delta_tildeN)
    attr(negLogL, "hessian") <- - hesMat
    ## return
    negLogL
}


## compute baseline hazard rate function
ic_h0t <- function(dat, tied) {
    numer <- ic_deltaTildeN(dat, tied)
    denom <- ic_k0(dat, tied)
    ifelse(denom > 0, numer / denom, 0)
}


## building blocks that follows notation in manuscript
ic_deltaTildeN <- function(dat, tied) {
    out <- with(dat, as.numeric(eventIdx) * p_jk)
    if (tied)
        return(aggregateSum(out, dat$time))
    out
}


## baseline hazard rate for censoring time
ic_h_c <- function(dat, tied) {
    numer <- ic_deltaC(dat, tied)
    denom <- revcumsum(dat$p_jk)
    if (tied)
        denom <- denom[dat$firstIdx]
    ifelse(numer > 0, numer / denom, 0)
}


ic_deltaC <- function(dat, tied) {
    out <- with(dat, as.numeric(! eventIdx) * p_jk)
    if (tied)
        return(aggregateSum(out, dat$time))
    out
}


ic_k0 <- function(dat, tied) {
    p_jk_xExp <- with(dat, p_jk * xExp)
    out <- revcumsum(p_jk_xExp)
    if (tied)
        out <- out[dat$firstIdx]
    out
}


ic_k1 <- function(parSeq, dat, xMat, tied) {
    ## matrix of dimension # unique event time by # parameters
    out <- sapply(parSeq, function(ind, dat, xMat) {
        p_jk_xExp_x <- with(dat, p_jk * xExp) * xMat[, ind]
        revcumsum(p_jk_xExp_x)
    }, dat = dat, xMat = xMat)
    if (tied)
        out <- out[dat$firstIdx, ]
    out
}


ic_k2 <- function(parSeq, dat, xMat, tied) {
    ind_grid <- expand.grid(parSeq, parSeq)
    ## matrix of dimension # unique event time by (# parameters) ^ 2
    out <- mapply(function(ind1, ind2) {
        p_jk_xExp_x2 <- with(dat, p_jk * xExp) * xMat[, ind1] * xMat[, ind2]
        revcumsum(p_jk_xExp_x2)
    }, ind_grid[, 1L], ind_grid[, 2L])
    if (tied)
        out <- out[dat$firstIdx, ]
    out
}


ic_dLbeta <- function(xMatDeltaN, k_0, k_1, delta_tildeN) {
    sum_jk <- colSums(xMatDeltaN)
    int_t <- colSums(k_1 / k_0 * delta_tildeN, na.rm = TRUE)
    sum_jk - int_t
}


ic_d2Lbeta <- function(parSeq, k_0, k_1, k_2, delta_tildeN) {
    ## part 1
    nPar <- length(parSeq)
    mat1 <- k_2 / k_0 * delta_tildeN
    part1 <- matrix(colSums(mat1, na.rm = TRUE), nPar, nPar)

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
ic_initPi <- function(censorRate, dat, equally = FALSE, ...)
{
    ## nonsense to eliminate cran checking note
    eventIdx <- NULL

    ## mixture probability for each subject: piVec
    numTab <- table(dat$ID)
    dat$numRecord <- numTab[match(as.character(dat$ID), names(numTab))]

    ## subject ID with censoring records
    cenID <- with(base::subset(dat, ! eventIdx), unique(ID))
    cenIdx <- as.integer(dat$ID %in% cenID)

    ## for subjects with multiple records
    idx1 <- with(dat, numRecord > 1L & eventIdx)
    idx2 <- with(dat, numRecord > 1L & (! eventIdx))

    piVec <- rep(1, NROW(dat))
    if (equally) {
        piVec[idx1 | idx2] <- 1 / dat[idx1 | idx2, "numRecord"]
    } else {
        piVec[idx1] <- (1 - censorRate * cenIdx[idx1]) /
            (dat[idx1, "numRecord"] - cenIdx[idx1])
        piVec[idx2] <- censorRate
    }
    piVec
}


## assign initial piVec not based on censoring
ic_initPi2 <- function(pi0, dat, randomly = FALSE, ...)
{
    nData <- NROW(dat)
    piVec <- rep(1L, nData)
    dupIdx <- duplicated(dat$ID)
    dupID <- dat$ID[dupIdx]
    uniIdx <- ! dat$ID %in% unique(dupID)
    numTab <- table(dupID)
    numRecord <- numTab[match(as.character(dupID), names(numTab))]
    firstIdx <- ! (dupIdx | uniIdx)
    piVec[firstIdx] <- pi0
    piVec[dupIdx] <- (1 - pi0) / numRecord
    ## randomly impute within subject
    if (randomly) {
        randomOrder <- stats::runif(nData)
        piVec <- piVec[order(dat$ID, randomOrder)]
    }
    ## return
    piVec
}


## \item \code{censorRate}: A positive numeric vector taking values between 0
##     and 1. If a vector of length more than one is given, multiple starting
##     values on prior probabilities are applied. In each initialization,
##     individual value in the vector is taken as parameter in the multinomial
##     distribution indicating the prior probability of censoring record being
##     true.  The remaining probability mass is equally assigned to the
##     possible event records. For example, when \code{censorRate = 1} is
##     specified, the extended model reduces to the regular Cox model fitting
##     data before matching, which takes all the censoring records as true. On
##     the contrary, when \code{censorRate = 0} is specified, the model does
##     not take possible censoring records into account.  If \code{multiStart =
##     FALSE}, the default value of \code{censorRate} is the sample censoring
##     rate among subjects having records without uncertainty.
##     If \code{parametricOnly}
##
## \item \code{multiStart}: A logical value specifying whether multiple
##     starting values of \code{censorRate} is needed. The argument is ignored
##     when explicit \code{censorRate} or \code{piVec} is specified. The
##     default value is \code{FALSE} for less computation time. However,
##     enabling multiple starting values \code{multiStart = TRUE} is suggested
##     to reduce the influence of possible local maximizer. If \code{TRUE},
##     \code{censorRate} takes a grid from 0 to 1 with step 0.02 by default.
##
## }
##
## take care of start values internally
iCoxph_start <- function(beta_vec = NULL,
                         beta_mat = NULL,
                         methods = NULL,
                         ...,
                         nBeta_,
                         dat_)
{
    ## nonsense to eliminate cran checking note
    ID <- NULL

    dupID <- with(dat_, unique(ID[duplicated(ID)]))
    dupIdx <- dat_$ID %in% dupID
    uniDat <- base::subset(dat_, ! dupIdx)
    uni_xMat <- as.matrix(uniDat[, - seq_len(3L)])
    xMat <- as.matrix(dat_[, - seq_len(3L)])

    ## censorRate0 <- round(1 - mean(uniDat$event), 2)
    ## if (is.null(censorRate)) {
    ##     ## set censorRate from sample truth data
    ##     ## if missing at random, the true censoring rate
    ##     ## can be estimated by true data of unique records
    ##     step_by <- 0.02
    ##     censorRate <- if (multiStart) {
    ##                       seq.int(0, 1, step_by)
    ##                   } else if (semiparametric_only) {
    ##                       NA_real_
    ##                   } else {
    ##                       censorRate0
    ##                   }
    ## } else if (any(censorRate > 1 | censorRate < 0)) {
    ##     stop("Starting prob. of censoring case being true",
    ##          "should between 0 and 1.", call. = FALSE)
    ## }

    ## use parametric estimate as starting values
    ## if (any(parametric)) {
    ##     event_funs <- ic_parametric_start(uniDat$time, uniDat$event, uni_xMat)
    ##     cen_funs <- ic_parametric_start(uniDat$time, 1 - uniDat$event)
    ##     hVec <- event_funs$haz_fun(dat_$time, xMat)
    ##     sVec <- event_funs$surv_fun(dat_$time, xMat)
    ##     h_cVec <- cen_funs$haz_fun(dat_$time)
    ##     G_cVec <- cen_funs$surv_fun(dat_$time)
    ##     if (parametric[1L]) {
    ##         ## following the equations derived
    ##         log_w_jk_1 <- ifelse(
    ##             dat_$event > 0,
    ##             log(hVec) + log(sVec) + log(G_cVec),
    ##             log(h_cVec) + log(G_cVec) + log(sVec)
    ##         )
    ##         w_jk_1 <- exp(log_w_jk_1)
    ##         p_jk_denom_1 <- aggregateSum(w_jk_1, dat_$ID, simplify = FALSE)
    ##         pi_par_1 <- ifelse(dupIdx, w_jk_1 / p_jk_denom_1, 1)
    ##         piMat <- cbind(piMat, pi_par_1)
    ##     }
    ##     if (length(parametric) > 1L && parametric[2L]) {
    ##         ## if hazard estimates is not reliable
    ##         log_w_jk_2 <- log(sVec) + log(G_cVec)
    ##         w_jk_2 <- exp(log_w_jk_2)
    ##         p_jk_denom_2 <- aggregateSum(w_jk_2, dat_$ID, simplify = FALSE)
    ##         pi_par_2 <- ifelse(dupIdx, w_jk_2 / p_jk_denom_2, 1)
    ##         piMat <- cbind(piMat, pi_par_2)
    ##     }
    ## }

    ## use non/semi-parametric estimates
    if (is.null(methods)) {
        stop("No iniliaization method was specified.", call. = FALSE)
    }
    pi_mat <- NULL
    ## initialization of baseline functions
    event_funs <- ic_semi_parametric_start(
        uniDat$time, uniDat$event, uni_xMat
    )
    cen_funs <- ic_semi_parametric_start(uniDat$time, 1 - uniDat$event)
    hVec <- event_funs$haz_fun(dat_$time, xMat)
    sVec <- event_funs$surv_fun(dat_$time, xMat)
    h_cVec <- cen_funs$haz_fun(dat_$time)
    G_cVec <- cen_funs$surv_fun(dat_$time)
    ## initialized baseline hazards with nearest (left) neighbor
    if ("nearest_hazard" %in% methods) {
        ## following the equations derived
        log_w_jk_3 <- ifelse(
            dat_$event > 0,
            log(hVec) + log(sVec) + log(G_cVec),
            log(h_cVec) + log(G_cVec) + log(sVec)
        )
        w_jk_3 <- exp(log_w_jk_3)
        p_jk_denom_3 <- aggregateSum(w_jk_3, dat_$ID, simplify = FALSE)
        pi_par_3 <- ifelse(dupIdx, w_jk_3 / p_jk_denom_3, 1)
        pi_mat <- cbind(pi_mat, "nearest_hazard" = pi_par_3)
    }
    ## use unit baseline hazards (or equivalently ignore them)
    if ("unit_hazard" %in% methods) {
        log_w_jk_4 <- log(sVec) + log(G_cVec)
        w_jk_4 <- exp(log_w_jk_4)
        p_jk_denom_4 <- aggregateSum(w_jk_4, dat_$ID, simplify = FALSE)
        pi_par_4 <- ifelse(dupIdx, w_jk_4 / p_jk_denom_4, 1)
        pi_mat <- cbind(pi_mat, "unit_hazard" = pi_par_4)
    }

    ## initialize covariate coefficient: beta
    if (is.null(beta_vec)) {
        ## if high censoring for subjects having unique records
        if (mean(uniDat$event) < 0.01) {
            beta_vec <- matrix(rep(0, nBeta_), ncol = 1)
        } else {
            uniDat$eventIdx <- NULL
            tmp <- tryCatch(
                with(uniDat, rcpp_coxph(time, event, uni_xMat)),
                warning = function(w) {
                    warning(w)
                    return(NULL)
                }, error = function(e) {
                    warning(e)
                    return(NULL)
                })
            beta <- if (is.null(tmp)) {
                        rep(0, nBeta_)
                    } else {
                        as.numeric(tmp$coef)
                    }
            beta_vec <- as.matrix(beta)
        }
    }
    beta_mat <- cbind(beta_mat, beta_vec)
    ## some quick checks on beta
    if (nrow(beta_mat) != nBeta_)
        stop(wrapMessages(
            "The number of starting values for coefficients of",
            "covariates does not match with the specified formula."
        ), call. = FALSE)
    ## some quick checks on pi
    if (! is.null(pi_mat)) {
        if (nrow(pi_mat) != nrow(dat_))
            stop("'pi_mat' must have same length with number of rows of data.")
        if (any(pi_mat > 1 | pi_mat < 0))
            stop("'pi_mat' has to be between 0 and 1.")
    }
    ## return
    list(beta_mat = beta_mat,
         pi_mat = pi_mat)
}


### return parametric hazard function and survival function
## ic_parametric_start <- function(time, event, xMat = NULL)
## {
##     ## transform estimates from AFT form to PH form for Weibull model
##     transCoef <- function(survRegObj) {
##         shape <- 1 / survRegObj$scale
##         betaEst <- survRegObj$coefficients
##         lambda0 <- exp(- shape * betaEst[1])
##         names(lambda0) <- NULL
##         betaEst <- - shape * betaEst[- 1]
##         list(beta = betaEst, lambda0 = lambda0)
##     }
##     ## fitting with survreg
##     fm <- if (is.null(xMat)) {
##               survival::Surv(time, event) ~ 1
##           } else {
##               survival::Surv(time, event) ~ xMat
##           }
##     fit <- survival::survreg(fm)
##     shape <- 1 / fit$scale
##     fit_list <- transCoef(fit)
##     betaHat <- fit_list$beta
##     lambda0 <- fit_list$lambda0
##     ## (baseline) hazard function and survival function
##     haz_fun <- function(time, xMat = NULL) {
##         betaX <- if (! is.null(xMat)) {
##                      as.numeric(xMat %*% betaHat)
##                  } else {
##                      0
##                  }
##         exp(log(lambda0) + log(shape) +
##             (shape - 1) * log(time) + betaX)
##     }
##     surv_fun <- function(time, xMat = NULL) {
##         betaX <- if (! is.null(xMat)) {
##                      as.numeric(xMat %*% betaHat)
##                  } else {
##                      0
##                  }
##         HVec <- exp(log(lambda0) + shape * log(time) + betaX)
##         exp(- HVec)
##     }
##     ## return
##     list(haz_fun = haz_fun,
##          surv_fun = surv_fun)
## }


### return semi-parametric hazard function and survival function
ic_semi_parametric_start <- function(time, event, xMat = NULL)
{
    ## set default value
    fit <- NULL
    if (! is.null(xMat)) {
        ## fitting may fail in bootstrap samples
        fit <- tryCatch(rcpp_coxph(time, event, xMat),
                        warning = function(w) {
                            return(NULL)
                        },
                        error = function(e) {
                            return(NULL)
                        })
    }
    if (is.null(fit)) {
        haz0 <- rcpp_mcf_right(time, event)
        haz0 <- data.frame(time = haz0$time,
                           hazard = haz0$cum_rate)
    } else {
        betaHat <- as.numeric(fit$coef)
        haz0 <- data.frame(time = fit$baseline$time,
                           hazard = fit$baseline$h0)
    }

    ## stepfun hazard function
    haz_fun <- function(time, xMat = NULL) {
        betaX <- if (! is.null(xMat)) {
                     as.numeric(xMat %*% betaHat)
                 } else {
                     0
                 }
        ## only consider positive hazard here
        inst_haz0 <- diff(c(0, haz0$hazard))
        inst_bool <- inst_haz0 > 0
        base_haz0 <- c(haz0$hazard[inst_bool][1], inst_haz0[inst_bool])
        base_haz_fun <- stats::stepfun(haz0$time[inst_bool], base_haz0)
        base_haz_fun(time) * exp(betaX)
    }

    ## survival function
    surv_fun <- function(time, xMat = NULL) {
        betaX <- if (! is.null(xMat)) {
                     as.numeric(xMat %*% betaHat)
                 } else {
                     0
                 }
        haz_fun <- stats::stepfun(haz0$time, c(0, haz0$hazard))
        exp(- haz_fun(time) * exp(betaX))
    }
    ## return
    list(haz_fun = haz_fun,
         surv_fun = surv_fun)
}
