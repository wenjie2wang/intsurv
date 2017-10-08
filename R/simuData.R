################################################################################
##
##   R package intsurv by Wenjie Wang, Kun Chen and Jun Yan
##   Copyright (C) 2017
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


##' Generate Simulated Data from Cox Model
##'
##' This function generates simulated survival data with uncertain records.
##' The baseline hazard function is parametric Weibull model
##' and the censoring times follow uniform distribution.
##'
##' @param nSubject Number of subjects
##' @param beta0Vec Time-invariant covariate coefficients.
##' @param xMat Design matrix. By default, three continuous variables following
##'     standard normal distribution and one binary variable following bernoulli
##'     distribution with equal probability are used.
##' @param maxNum Maximum number of uncertain records.
##' @param nRecordProb Probability of the number of uncertian records.
##' @param matchCensor The matching rate for subjects actually having censoring
##'     times.
##' @param matchEvent The matching rate for subjects actually having event
##'     times.
##' @param censorMin The lower boundary of the uniform distribution for
##'     generating censoring time.
##' @param censorMax The upper boundary of the uniform distribution for
##'     generating censoring time.
##' @param lambda A positive number, scale parameter in baseline rate function
##'     for true event times.
##' @param rho A positive number, shape parameter in baseline rate function for
##'     true event times.
##' @param fakeLambda1 A positive number, scale paramenter in baseline rate
##'     function for fake event times from one distribution.
##' @param fakeRho1 A positive number, shape parameter in baseline rate function
##'     for fake event times from one distribution.
##' @param fakeLambda2 A positive number, scale paramenter in baseline rate
##'     function for fake event times from another distribution.
##' @param fakeRho2 A positive number, shape parameter in baseline rate function
##'     for fake event times from another distribution.
##' @param mixture The mixture weights, i.e., the probabilies (summing up to
##'     one) of fake event times coming from different mixture components.
##' @param randomMiss A boolean value specifying whether the labels of the true
##'     records are missing completely at random (MCAR) or missing not at random
##'     (MNAR). The default value is \code{TRUE} for MCAR.
##' @param eventOnly A boolean value specifying whether the uncertain records
##'     only include possible events. The default value is \code{FALSE}, which
##'     considers the censoring cases as the possible truth in addition to event
##'     records.
##' @param ... Other arguments for future usage.
##'
##' @return A data frame with the following columns, \itemize{ \item \code{ID}:
##'     subject ID \item \code{obsTime}: observed event times \item
##'     \code{eventInd}: event indicators \item \code{latentInd}: latent label
##'     indicating the true recorods } and the corresponding covariates.
##'
##' @references
##'
##' Bender, R., Augustin, T. and Blettner, M. (2005),
##' Generating survival times to simulate Cox proportional hazards models.
##' Statistics in Medicine, 24: 1713--1723. doi:10.1002/sim.2059.
##'
##' Nadarajah, S., & Kotz, S. (2006). R programs for truncated
##' distributions. Journal of Statistical Software, 16(c02).
##'
##' @examples
##' dat <- simuWeibull(nSubject = 100)
##' @importFrom stats runif rnorm
##' @export
simuWeibull <- function(nSubject = 1e3, beta0Vec, xMat, maxNum = 2,
                        nRecordProb = c(0.9, 0.1),
                        matchCensor = 0.1, matchEvent = 0.1,
                        censorMin = 0.5, censorMax = 12.5,
                        lambda = 0.005, rho = 0.7,
                        fakeLambda1 = lambda * exp(- 3), fakeRho1 = rho,
                        fakeLambda2 = lambda * exp(3), fakeRho2 = rho,
                        mixture = 0.5, randomMiss = TRUE,
                        eventOnly = FALSE, ...)
{
    ## 1. Generate the true dataset with unique event time per subject and
    ## cencoring time. The censoring rate of the truth data set can be
    ## controled.

    ## 2. Take a sample of subject ID, where the sample size depends on the
    ## matching rate from original data to death database.

    ## 3. Generate fake records for selected the subjects.

    ## 3.1. If the subject has event time censored, i.e. the truth is the
    ## censoring record, generate given number of fake event times by different
    ## x's or distribution parameters and only keep fake events that happened
    ## before censoring time. (It is possible that not much fake events left,
    ## which means that the alternative distribution of fake events is probably
    ## not reasonable for the censoring time picked.)

    ## 3.2. If the subject has event time recorded, i.e. the truth is the event
    ## time, generate event times by different x's or distribution parameters.
    ## Take them and the censoring time picked as all the fake
    ## records. Similarly, only keep event times before the censoring time.

    ## non-sense just to surpress complaint from R check
    ID <- NULL

    ## set design matrix and coefficients for event times
    if (missing(xMat)) {
        x1 <- round(rnorm(nSubject), 3)
        x2 <- round(rnorm(nSubject), 3)
        x3 <- round(rnorm(nSubject), 3)
        x4 <- sample(c(0, 1), size = nSubject, replace = TRUE)
        xMat <- cbind(x1, x2, x3, x4)
    } else {
        xMat <- as.matrix(xMat)
    }
    beta0Vec <- if (missing(beta0Vec))
                    rep(1, ncol(xMat))
                else
                    as.numeric(beta0Vec)
    ## event times
    eventTime <- denWeib(nSubject, beta0Vec, xMat, lambda, rho, gen = TRUE)

    ## censoring times
    censorTime <- runif(n = nSubject, min = censorMin, max = censorMax)

    eventInd <- censorTime > eventTime
    obsTime <- ifelse(eventInd, eventTime, censorTime)
    trueDat <- data.frame(ID = seq_len(nSubject), xMat, obsTime = obsTime,
                          eventTime = eventTime, censorTime = censorTime,
                          eventInd, latentInd = 1L)

    ## match records for subjects actually having events
    numSubFake1 <- round(matchEvent * sum(trueDat$eventInd))
    subIDfake1 <- if (randomMiss) {
                      with(subset(trueDat, eventInd == 1L),
                           sample(ID, size = numSubFake1))
                  } else {
                      with(subset(trueDat, eventInd == 1L),
                           sample(ID, size = numSubFake1,
                                  prob = eventTime))
                  }

    ## match records for censoring cases
    numSubFake0 <- round(matchCensor * sum(! trueDat$eventInd))
    subIDfake0 <- if (randomMiss) {
                      with(subset(trueDat, eventInd == 0L),
                           sample(ID, size = numSubFake0))
                  } else {
                      with(subset(trueDat, eventInd == 0L),
                           sample(ID, size = numSubFake0,
                                  prob = eventTime))
                  }

    if (maxNum == 1) nRecordProb <- 1
    nRecords <- sample(maxNum, size = numSubFake0 + numSubFake1,
                       replace = TRUE, prob = nRecordProb)
    numFake <- sum(nRecords) - numSubFake1
    multiInd <- c(rep(c(subIDfake0, subIDfake1), nRecords - 1), subIDfake0)
    fakeDat <- trueDat[multiInd, ]
    xMat <- xMat[multiInd, ]

    ## fake event time from alternative dist. # 1 and # 2
    fakes1 <- denWeib(numFake, beta0Vec, xMat, lambda = fakeLambda1,
                      rho = fakeRho1, censorTime = fakeDat$censorTime,
                      gen = TRUE)
    fakes2 <- denWeib(numFake, beta0Vec, xMat, lambda = fakeLambda2,
                      rho = fakeRho2, censorTime = fakeDat$censorTime,
                      gen = TRUE)

    ## fake censored data for all selected subject with true event times
    fakeCensor <- subset(trueDat, ID %in% subIDfake1 & eventInd)
    if (eventOnly || ! nrow(fakeCensor)) {
        fakeCensor <- NULL
    } else {
        fakeCensor$obsTime <- fakeCensor$censorTime
        fakeCensor$eventInd <- fakeCensor$latentInd <- 0L
    }

    ## setup mixture of alternative density
    indDist1 <- runif(numFake)
    fakeDat$obsTime <- ifelse(indDist1 < mixture, fakes1, fakes2)
    fakeDat$eventInd <- 1L
    fakeDat$latentInd <- 0L

    ## remove all fake events after censorTime
    ## fakeDat <- subset(fakeDat, obsTime <= censorTime)
    ## by using truncated distribution, no need to remove

    ## combine all parts as output
    outDat <- rbind(trueDat, fakeDat, fakeCensor)
    outDat$eventInd <- as.integer(outDat$eventInd)
    out <- outDat[order(outDat$ID), ]
    row.names(out) <- out$eventTime <- out$censorTime <- NULL
    out
}


### internal function ==========================================================
denWeib <- function (nEvent = 1e3, beta0Vec, xMat, lambda, rho,
                     censorTime = NULL, gen = FALSE, ...) {

    if (missing(beta0Vec))
        beta0Vec <- c(1, 1)
    if (missing(xMat))
        xMat <- rep(0, length(beta0Vec))

    vVec <- runif(nEvent)
    if (gen) {
        expXbeta <- exp(xMat %*% beta0Vec)
        if (is.null(censorTime)) {
            eventTime <- (- log(1 - vVec) / (lambda * expXbeta)) ^ (1 / rho)
            return(eventTime)
        }
        ## generate event times before given censoring time
        ## using truncated distribution
        ## reference: Nadarajah and Kotz (2006). JSS.
        Gb <- 1 - exp(- lambda * censorTime ^ rho * expXbeta)
        eventTime <- (- log(1 - vVec * Gb) / (lambda * expXbeta)) ^ (1 / rho)
        idx <- eventTime >= censorTime
        return(eventTime)
    }

    ## otherwise
    newScale <- lambda * exp(xMat %*% beta0Vec)
    eventTime <- (- log(vVec) / newScale) ^ (1 / rho)
    tVec <- seq.int(0, max(eventTime), by = 0.01)
    htSt <- newScale * rho * tVec ^ (rho - 1) * exp(- newScale * tVec ^ rho)
    list(time = tVec, htSt = htSt, x = xMat, beta = beta0Vec)
}
