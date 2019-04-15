##' Simulating Data from Cox Cure Model with Uncertain Endpoints
##'
##' FIXME
##'
##' @examples
##' set.seed(808)
##' nSubject <- 1e2
##' coxMat <- cbind(rnorm(nSubject), rbinom(nSubject, size = 1, prob = 0.5))
##' cureMat <- cbind(rnorm(nSubject), rbinom(nSubject, size = 1, prob = 0.5), 1)
##' testDat <- simuCureUncertain(nSubject = nSubject, shape = 2, scale = 1e-3,
##'                              censorMin = 0, censorMax = 60,
##'                              coxCoef = c(0.5, 0.3), coxMat = coxMat,
##'                              cureCoef = c(0.3, 0.8, - 1), cureMat = cureMat)
##'
simuCureUncertain <- function(nSubject = 1e3, shape = 2, scale = 0.001,
                              censorMin = 0, censorMax = 60,
                              coxCoef, coxMat,
                              cureCoef, cureMat = coxMat,
                              n1_rate = 0.5,
                              n21_rate = 0.5,
                              n22_rate = 0.5,
                              n31_rate = 1 - n1_rate,
                              n32_rate = 1 - n21_rate,
                              n33_rate = 1 - n22_rate,
                              type = c("weibull", "gompertz"),
                              ...)
{
    ## a bare-bones version without any checking and pre-processing on inputs

    ## 0. generate subjects' ID
    id <- seq_len(nSubject)

    ## 1. generate cure indicator for each subject based on logistics model.
    cure_p <- 1 / (1 + exp(- as.numeric(cureMat %*% cureCoef)))
    cure_ind <- rbinom(nSubject, size = 1, prob = cure_p) < 1
    num_risk_subjects <- sum(! cure_ind)

    ## 2. generate event times, censoring times independently, and thus observed
    ##    times and event indicator for subjects that are not cured.
    expXbeta <- as.numeric(exp(coxMat[! cure_ind, ] %*% coxCoef))
    type <- match.arg(type)
    eventTimes <- if (type == "gompertz") {
                      rGompertz(num_risk_subjects, shape, scale * expXbeta)
                  } else {
                      rWeibull(num_risk_subjects, shape, scale * expXbeta)
                  }
    censorTimes <- runif(num_risk_subjects, min = censorMin, max = censorMax)
    obsTimes <- pmin(eventTimes, censorTimes)
    event_ind <- eventTimes < censorTimes

    ## 3. for those subjects (n_{22} + n_{33}) that are actually cured, generate
    ##    censoring times for them and randomly assign n_{22} of them to case
    ##    2b, and n_{33} of them to case 3c.
    id_cure <- id[cure_ind]
    num_cure <- nSubject - num_risk_subjects
    censorTimesCure <- runif(num_cure, min = censorMin, max = censorMax)
    n22 <- floor(num_cure * n22_rate)
    n33 <- nSubject - num_risk_subjects - n22
    id_n22 <- sample(id_cure, size = n22)
    id_n33 <- setdiff(id_cure, id_n22)
    n22_dat <- if (n22 > 0) {
                   data.frame(
                       obs_time = censorTimesCure[seq_len(n22)],
                       obs_event = 0, oracle_event = 0,
                       oracle_cure = 1, case = "n_22", id = id_n22
                   )
               } else {
                   NULL
               }
    n33_dat <- if (n33 > 0) {
                   data.frame(
                       obs_time = censorTimesCure[- seq_len(n22)],
                       obs_event = NA_integer_, oracle_event = 0,
                       oracle_cure = 1, case = "n_33", id = id_n33
                   )
               } else {
                   NULL
               }

    ## 4. sample n_{31} of (n_{1} + n_{31}) subjects that are actually having
    ##    events and add fake censoring times for them, where the fake censoring
    ##    times are simply the true event times.
    id_event <- id[! cure_ind][event_ind]
    num_event <- sum(event_ind)
    n1 <- floor(num_event * n1_rate)
    n31 <- num_event - n1
    id_n1 <- sort(sample(id_event, size = n1))
    id_n31 <- sort(setdiff(id_event, id_n1))
    n1_dat <- if (n1 > 0) {
                  data.frame(
                      obs_time = obsTimes[event_ind][id_event %in% id_n1],
                      obs_event = 1, oracle_event = 1,
                      oracle_cure = 0, case = "n_1", id = id_n1
                  )
              } else {
                  NULL
              }
    n31_dat <- if (n31 > 0) {
                   data.frame(
                       obs_time = obsTimes[event_ind][id_event %in% id_n31],
                       obs_event = NA_integer_, oracle_event = 1,
                       oracle_cure = 0, case = "n_31", id = id_n31
                   )
               } else {
                   NULL
               }

    ## 5. sample n_{32} of (n_{21} + n_{32}) subjects that are actually censored
    ##    and add fake event times for them, where the fake event times are
    ##    simply the true censoring times.
    id_censor <- id[! cure_ind][! event_ind]
    num_censor <- num_risk_subjects - num_event
    n21 <- floor(num_censor * n21_rate)
    n32 <- num_censor - n21
    id_n21 <- sort(sample(id_censor, size = n21))
    id_n32 <- sort(setdiff(id_censor, id_n21))
    n21_dat <- if (n21 > 0) {
                   data.frame(
                       obs_time = obsTimes[! event_ind][id_censor %in% id_n21],
                       obs_event = 0, oracle_event = 0,
                       oracle_cure = 0, case = "n_21", id = id_n21
                   )
               } else {
                   NULL
               }
    n32_dat <- if (n32 > 0) {
                   data.frame(
                       obs_time = obsTimes[! event_ind][id_censor %in% id_n32],
                       obs_event = NA_integer_, oracle_event = 0,
                       oracle_cure = 0, case = "n_32", id = id_n32
                   )
               } else {
                   NULL
               }

    ## combine all the simulated data together
    outDat <- rbind(n1_dat, n21_dat, n22_dat, n31_dat, n32_dat, n33_dat)
    outDat <- outDat[order(outDat$id), ]
    xDat <- data.frame(coxMat)
    colnames(xDat) <- paste0("x", seq_len(ncol(xDat)))
    zDat <- data.frame(cureMat)
    colnames(zDat) <- paste0("z", seq_len(ncol(zDat)))
    xzDat <- cbind(id = id, xDat, zDat)
    outDat <- merge(outDat, xzDat, by = "id")

    ## add some attributes
    attr(outDat, "n") <- setNames(
        c(nSubject, n1, n21 + n22, n21, n22, n31 + n32 + n33, n31, n32, n33),
        c("n", "n1", "n2", "n21", "n22", "n3", "n31", "n32", "n33")
    )

    ## return
    outDat
}


##' Generate Simulated Survival Data
##'
##' The function generates survival data with exact event times and
##' right-censoring time.  The event times can be simulated from Gompertz
##' distribution, Weibull distribution, and exponential distribution as a
##' special case of Weibull distribution (\code{shape = 1}).  Censoring times
##' are generated from uniform distribution.
##'
##' @param nSubject Number of subjects.
##' @param beta0 Time-invariant covariate coefficients.
##' @param xMat Design matrix.
##' @param censorMin A non-negative number, lower bound of uniform
##'     distribution for censoring times.
##' @param censorMax A positive number, upper bound of uniform distribution
##'     for censoring times.
##' @param shape A positive number, shape parameter in baseline hazard
##'     function for event times.
##' @param scale A positive number, scale parameter in baseline hazard
##'     function for event times.
##' @param type Distribution name of event times. Partial matching is allowed.
##' @return A data frame.
##' @author Wenjie Wang
##' @references
##' Bender, R., Augustin, T., & Blettner, M. (2005).
##' Generating survival times to simulate Cox proportional hazards models.
##' \emph{Statistics in Medicine}, 24(11), 1713--1723.
##' @examples
##' simuData(100, type = "weibull")
simuData <- function(nSubject = 1e3, beta0, xMat, censorMin = 0,
                     censorMax = 10, shape = 1, scale = 0.01,
                     type = c("gompertz", "weibull"), ...) {

    ## set design matrix and coefficients for event times if missing
    if (missing(xMat)) {
        x1 <- sample(c(0, 1), size = nSubject, replace = TRUE)
        x2 <- round(rnorm(nSubject), 3)
        xMat <- cbind(x1, x2)
    } else
        xMat <- as.matrix(xMat)

    beta0 <- if (missing(beta0)) c(1, 1) else as.numeric(beta0)
    expXbeta <- as.numeric(exp(xMat %*% beta0))

    ## censoring times
    censorTime <- runif(n = nSubject, min = censorMin, max = censorMax)

    ## event times
    type <- match.arg(type)
    eventTime <- if (type == "gompertz")
                     rGompertz(nSubject, shape, scale * expXbeta)
                 else
                     rWeibull(nSubject, shape, scale * expXbeta)

    ## event indicator
    status <- as.integer(eventTime < censorTime)

    ## observed times
    obsTime <- ifelse(status, eventTime, censorTime)

    ## return
    data.frame(ID = seq_len(nSubject), xMat, obsTime, status)
}


## simulate one by one
simuIntCure <- function(nSubject = 1e3, shape = 2, scale = 0.001,
                        censorMin = 0, censorMax = 60,
                        p1 = 0.5, p2 = 0.5, p3 = 0.5,
                        coxCoef, coxMat,
                        cureCoef, cureMat = coxMat,
                        type = c("weibull", "gompertz"),
                        ...)
{
    type <- match.arg(type)

    simuIntCureOne <- function(i)
    {
        ## 1. generate cure indicator for each subject based on logistics model.
        cure_p <- 1 / (1 + exp(- as.numeric(cureMat[i, ] %*% cureCoef)))
        cure_ind <- rbinom(1, size = 1, prob = cure_p) < 1

        ## 2.1 if cured
        if (cure_ind) {
            censorTime <- runif(1, min = censorMin, max = censorMax)
            b2 <- rbinom(1, 1, p2)
            oracle_case <- ifelse(b2, "2_b", "3_c")
            obs_event <- ifelse(b2, 0, NA)
            out <- data.frame(obs_time = censorTime, obs_event = NA,
                              oracle_event = 0, oracle_cure = 1,
                              case = "3_c")
            return(out)
        }

        ## 2.2 otherwise (not cured)
        expXbeta <- as.numeric(exp(coxMat[i, ] %*% coxCoef))
        eventTime <- if (type == "gompertz") {
                         rGompertz(1, shape, scale * expXbeta)
                     } else {
                         rWeibull(1, shape, scale * expXbeta)
                     }
        censorTime <- runif(1, min = censorMin, max = censorMax)
        obsTime <- min(eventTime, censorTime)
        obsEvent <- as.integer(eventTime < censorTime)
        if (obsEvent) {
            ## case 1 and case 3a
            b1 <- rbinom(1, 1, p1)
            oracle_case = ifelse(b1, "1", "3_a")
            obs_event <- ifelse(b1, obsEvent, NA)
        } else {
            ## case 2_a and case 3b
            b3 <- rbinom(1, 1, p3)
            oracle_case = ifelse(b3, "2_a", "3_b")
            obs_event <- ifelse(b3, obsEvent, NA)
        }
        out <- data.frame(obs_time = obsTime, obs_event = obsEvent,
                          oracle_event = obsEvent, oracle_cure = 0,
                          case = oracle_case)
        return(out)
    }
    res <- do.call(rbind, lapply(seq_len(nSubject), simuIntCureOne))
    colnames(coxMat) <- paste0("x", seq_len(ncol(coxMat)))
    colnames(cureMat) <- paste0("z", seq_len(ncol(cureMat)))
    cbind(res, coxMat, cureMat)
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
