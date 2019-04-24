Rcpp::sourceCpp("../src/int_coxph_cure.cpp")

source("simu-data.R")

set.seed(808)

nSubject <- 1e4
coxMat <- cbind(rnorm(nSubject), rnorm(nSubject))
cureMat <- cbind(1, rnorm(nSubject), rnorm(nSubject))

testDat <- simuIntCure_all(nSubject_case12 = nSubject * 0.8,
                           nSubject_case3 = nSubject * 0.2,
                           event_shape = 2, event_scale = 0.1,
                           censor_shape = 2.5, censor_scale = 0.001,
                           prob_event = 0.5,
                           coxCoef = c(0.5, 0.3), coxMat = coxMat,
                           cureCoef = c(0, 2, - 2), cureMat = cureMat)

xx <- sort(testDat$obs_time)
plot(xx, oracle_d0(xx, 2, 0.1), type = "l", col = "red")
lines(xx, oracle_d0(xx, 2.5, 0.001), type = "l")

testDat$case <- factor(testDat$case,
                       levels = sort(unique(as.character(testDat$case))))

table(testDat$case)
xtabs(~ obs_event + case, data = testDat, addNA = TRUE)

testDat$obs_event2 <- testDat$obs_event
testDat$obs_event2[is.na(testDat$obs_event2)] <- 0.01

## sort the testDat by time
testDat <- testDat[order(testDat$obs_time, 1 - testDat$obs_event2), ]

testDat$obs_time2 <- round(testDat$obs_time)
## testDat$obs_event[is.na(testDat$obs_event)] <- 1

table(testDat$obs_time2)
table(testDat$obs_time2[testDat$obs_event == 0])
table(testDat$obs_time2[testDat$obs_event == 1])
table(testDat$obs_time2[is.na(testDat$obs_event)])

tmp <- with(testDat,
            int_coxph_cure(obs_time, obs_event2,
                           cbind(x1, x2), cbind(z2, z3),
                           cox_mstep_rel_tol = 1e-4,
                           cox_mstep_max_iter = 100,
                           cure_mstep_rel_tol = 1e-4,
                           cure_mstep_max_iter = 100,
                           prob_event_start = 0.5,
                           cox_start = c(0.5, 0.3),
                           cure_start = c(0, 2, -2),
                           em_max_iter = 200,
                           em_rel_tol = 1e-4))

tmp$cure_beta
tmp$cox_beta

boxplot(tmp$w_mat)

hist(tmp$prob_not_cured_case3)
## plot(xx, tmp$h0 * tmp$S0, col = "red")
## points(xx, tmp$hc * tmp$Sc)


## posterior probs
prob_cured <- as.numeric(1 - tmp$prob_not_cured_case3)
prob_event_not_cured <- tmp$prob_event * (1 - prob_cured)
prob_censored_not_cured <- (1 - tmp$prob_event) * (1 - prob_cured)
## summary(prob_cured + prob_event_not_cured + prob_censored_not_cured)
case3 <- testDat$case[substr(testDat$case, 1, 1) == "3"]

case_3a_ind <- which(case3 == "3_a")
summary(prob_event_not_cured[case_3a_ind])
hist(prob_event_not_cured[case_3a_ind], breaks = 30)

case_3b_ind <- which(case3 == "3_b")
summary(prob_censored_not_cured[case_3b_ind])
hist(prob_censored_not_cured[case_3b_ind], breaks = 30)

y_3ab <- 2 - as.numeric(factor(testDat$case[testDat$case %in% c("3_a", "3_b")]))
x_3ab <- prob_event_not_cured[case3 %in% c("3_a", "3_b")]
## pROC::roc(y_3ab, x_3ab)

case_3c_ind <- which(case3 == "3_c")
hist(prob_cured[case_3c_ind], breaks = 30)

case_2a_ind <- which(testDat$case == "2_a")
hist(tmp$prob_not_cured[case_2a_ind], breaks = 30)

case_2b_ind <- which(testDat$case == "2_b")
hist(1 - tmp$prob_not_cured[case_2b_ind], breaks = 30)


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
