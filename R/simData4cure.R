##
## intsurv: Integrative Survival Models
## Copyright (C) 2017-2020  Wenjie Wang <wjwang.stat@gmail.com>
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


##' Simulate Data from Cox Cure Model with Uncertain Event Status
##'
##' @param nSubject A positive integer specifying number of subjects.
##' @param shape A positive number specifying the shape parameter of the
##'     distribution of the event times.
##' @param scale A positive number specifying the scale parameter of the
##'     distribution of the event times.
##' @param lambda_censor A positive number specifying the rate parameter of the
##'     exponential distribution for generating censoring times.
##' @param max_censor A positive number specifying the largest censoring time.
##' @param p1 A number between 0 and 1 specifying the probability of simulating
##'     events with observed event indicators given the simulated event times.
##' @param p2 A number between 0 and 1 specifying the probability of simulating
##'     susceptible censoring times with observed event status given the
##'     simulated susceptible censoring times.
##' @param p3 A number between 0 and 1 specifying the probability of simulating
##'     cured censoring times with observed event status given the simulated
##'     cured censoring times.
##' @param survMat A numeric matrix representing the design matrix of the
##'     survival model part.
##' @param cureMat A numeric matrix representing the design matrix excluding
##'     intercept of the cure rate model part.
##' @param b0 A number representing the intercept term for the cure rate model
##'     part.
##' @param survCoef A numeric vector for the covariate coefficients of the
##'     survival model part.
##' @param cureCoef A numeric vector for the covariate coefficients of the
##'     cure model part.
##' @param ... Other arguments not used currently.
##'
##' @return
##' A data.frame with the following columns:
##' \itemize{
##'
##' \item \code{obs_time}: Observed event/survival times.
##'
##' \item \code{obs_event}: Observed event status.
##'
##' \item \code{event_time}: Underlying true event times.
##'
##' \item \code{censor_time}: underlying true censoring times.
##'
##' \item \code{oracle_event}: underlying true event indicators.
##'
##' \item \code{oracle_cure}: underlying true cure indicators.
##'
##' \item \code{case}: underlying true case labels.
##'
##' }
##'
##' @references
##'
##' Wang, W., Luo, C., Aseltine, R. H., Wang, F., Yan, J., & Chen, K. (2020).
##' Suicide Risk Modeling with Uncertain Diagnostic Records. \emph{arXiv
##' preprint arXiv:2009.02597}.
##'
##' @examples
##' ## see examples of function cox_cure
##' @importFrom stats binomial rbinom rexp runif
##' @export
simData4cure <- function(nSubject = 1e3,
                         shape = 2, scale = 0.1,
                         lambda_censor = 1, max_censor = Inf,
                         p1 = 0.9, p2 = 0.9, p3 = 0.9,
                         survMat, cureMat = survMat,
                         b0 = stats::binomial()$linkfun(0.7),
                         survCoef = rep(1, ncol(survMat)),
                         cureCoef = rep(1, ncol(cureMat)),
                         ...)
{
    simu_one <- function(i)
    {
        ## 1. generate cure indicator for each subject based on logistics model.
        expit <- binomial()$linkinv
        not_cure_p <- expit(as.numeric(cureMat[i, ] %*% cureCoef + b0))
        cure_ind <- rbinom(1, size = 1, prob = not_cure_p) < 1
        expXbeta <- as.numeric(exp(survMat[i, ] %*% survCoef))

        if (cure_ind) {
            ## 2.1 if cured
            censorTime <- min(rexp(1, rate = lambda_censor), max_censor)
            b <- rbinom(1, 1, p3)
            out <- data.frame(obs_time = censorTime,
                              obs_event = ifelse(b, 0, NA),
                              event_time = Inf,
                              censor_time = censorTime,
                              oracle_event = 0,
                              oracle_cure = 1,
                              case = ifelse(b, "2b", "3c"))
        } else {
            ## 2.2 otherwise (not cured)
            eventTime <- rWeibull(1, shape, scale * expXbeta)
            censorTime <- min(rexp(1, rate = lambda_censor), max_censor)
            obsTime <- min(eventTime, censorTime)
            obsEvent <- as.integer(eventTime <= censorTime)
            if (obsEvent) {
                b <- rbinom(1, 1, p1)
                out <- data.frame(obs_time = obsTime,
                                  obs_event = ifelse(b, 1, NA),
                                  event_time = eventTime,
                                  censor_time = censorTime,
                                  oracle_event = 1,
                                  oracle_cure = 0,
                                  case = ifelse(b, "1", "3a"))
            } else {
                b <- rbinom(1, 1, p2)
                out <- data.frame(obs_time = obsTime,
                                  obs_event = ifelse(b, 0, NA),
                                  event_time = eventTime,
                                  censor_time = censorTime,
                                  oracle_event = 0,
                                  oracle_cure = 0,
                                  case = ifelse(b, "2a", "3b"))
            }
        }
        ## Bayes classification
        ## event_den <- dWeibull(out$obs_time, shape, scale * expXbeta)
        ## censor_den <- dexp(out$obs_time, lambda_censor)
        ## event_surv <- 1 - pWeibull(out$obs_time, shape, scale * expXbeta)
        ## censor_surv <- 1 - pexp(out$obs_time, lambda_censor)
        ## m1 <- event_den * censor_surv * not_cure_p
        ## m2 <- censor_den * event_surv * not_cure_p
        ## m3 <- (1 - not_cure_p) * censor_den
        ## out$bayes_event <- m1 / (m1 + m2 + m3)
        ## return
        return(out)
    }
    res <- do.call(rbind, lapply(seq_len(nSubject), simu_one))
    colnames(survMat) <- paste0("x", seq_len(ncol(survMat)))
    colnames(cureMat) <- paste0("z", seq_len(ncol(cureMat)))
    cbind(res, survMat, cureMat)
}


### internal function ==========================================================
## similar function with eha::rgompertz, where param == "canonical"
rGompertz <- function(n, shape = 1, scale = 1, ...) {
    u <- runif(n)
    1 / shape * log1p(- shape * log1p(- u) / (scale * shape))
}
## similar function with stats::rweibull but with different parametrization
## reduces to exponential distribution when shape == 1
rWeibull <- function(n, shape = 1, scale = 1, ...) {
    u <- runif(n)
    (- log1p(- u) / scale) ^ (1 / shape)
}
hWeibull <- function(x, shape = 1, scale = 1) {
    scale * shape * x ^ (shape - 1)
}
HWeibull <- function(x, shape = 1, scale = 1) {
    scale * x ^ shape
}
pWeibull <- function(x, shape = 1, scale = 1) {
    1 - exp(- scale * x ^ shape)
}
dWeibull <- function(x, shape = 1, scale = 1) {
    hWeibull(x, shape, scale) * exp(- scale * x ^ shape)
}
