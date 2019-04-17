Rcpp::sourceCpp("../src/int_coxph_cure.cpp")

source("simu-data.R")

set.seed(808)

nSubject <- 300
coxMat <- cbind(rnorm(nSubject), rnorm(nSubject))
cureMat <- cbind(1, rnorm(nSubject), rnorm(nSubject))

testDat <- simuIntCure(nSubject = nSubject, shape = 1, scale = 0.4,
                       censorMin = 1, censorMax = 10, event_pi = 0.5,
                       p1 = 0.5, p2 = 0.5, p3 = 0.5,
                       coxCoef = c(0.5, 0.3), coxMat = coxMat,
                       cureCoef = c(0, 0.5, - 0.5), cureMat = cureMat)

testDat$case <- factor(testDat$case,
                       levels = sort(unique(as.character(testDat$case))))
table(testDat$case)

testDat$obs_event2 <- testDat$obs_event
testDat$obs_event2[is.na(testDat$obs_event2)] <- 0.5

testDat$obs_time2 <- round(testDat$obs_time)
## testDat$obs_event[is.na(testDat$obs_event)] <- 1

table(testDat$obs_time2)
table(testDat$obs_time2[testDat$obs_event == 0])
table(testDat$obs_time2[testDat$obs_event == 1])

tmp <- with(testDat,
            int_coxph_cure(obs_time, obs_event2,
                           cbind(x1, x2), cbind(z2, z3),
                           cox_mstep_rel_tol = 1e-4,
                           cox_mstep_max_iter = 100,
                           cure_mstep_rel_tol = 1e-4,
                           cure_mstep_max_iter = 100,
                           prob_event_start = 0.5,
                           cox_start = c(0, 0),
                           cure_start = c(0.1, 1, - 1),
                           em_max_iter = 1000,
                           em_rel_tol = 1e-4))


## 2. simulate some right censoring data
library(reda)
set.seed(123)
n <- 100
p <- 3
z2 <- matrix(rnorm(n * p), nrow = n, ncol = p)
testDat <- simEventData(z = z2, endTime = 10,
                        rho = 0.1, zCoef = c(0.1, 0.2, - 0.1),
                        recurrent = FALSE)
testDat$time2 <- round(testDat$time)
with(testDat, table(event, time2))
table(testDat$event)

idx <- sample(which(testDat$event > 0), size = 30)
testDat$event[idx] <- 0.5
table(testDat$event)
with(testDat, table(event, time2))

tmp <- with(testDat,
            int_coxph_cure(time2, event, z2, z2,
                           cox_mstep_rel_tol = 1e-6,
                           cox_mstep_max_iter = 100,
                           cure_mstep_rel_tol = 1e-6,
                           cure_mstep_max_iter = 100,
                           prob_event_start = 0.5,
                           cox_start = c(0, 0),
                           cure_start = c(0, 0, 0, 0),
                           em_max_iter = 30,
                           em_rel_tol = 1e-4))


## stepwise debug
