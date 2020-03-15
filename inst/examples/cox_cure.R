library(intsurv)

### regular Cox cure rate model ======================================
## 1. simulate right-censored data with a cure fraction
set.seed(123)
n_obs <- 2e2
p <- 5
x_mat <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
colnames(x_mat) <- paste0("x", seq_len(p))
cure_beta <- rep(0.5, p)
b0 <- - 1
expit <- binomial()$linkinv
ncure_prob <- expit(as.numeric(b0 + x_mat %*% cure_beta))
is_cure <- 1 - rbinom(n_obs, size = 1, prob = ncure_prob)
surv_beta <- rep(0.5, p)
risk_score <- as.numeric(x_mat %*% surv_beta)
event_time <- rexp(n_obs, exp(as.numeric(x_mat %*% surv_beta)))
censor_time <- 10
event <- ifelse(event_time < censor_time & ! is_cure, 1, 0)
obs_time <- ifelse(event > 0, event_time, censor_time)

## model-fitting from given design matrices
fit1 <- cox_cure.fit(x_mat, x_mat, obs_time, event, bootstrap = 30)
summary(fit1)

## coefficient estimates from both model parts
coef(fit1)

## or a particular part
coef(fit1, "surv")
coef(fit1, "cure")

## weighted concordance index (C-index)
fit1$model$c_index
## which also can be computed as follows
cIndex(time = obs_time, event = event,
       risk_score = fit1$fitted$surv_xBeta,
       weight = ifelse(event > 0, 1, fit1$fitted$susceptible_prob))

## 2. create a toy dataset
toy_dat <- data.frame(time = obs_time, status = event)
toy_dat$group <- cut(abs(x_mat[, 1L]), breaks = c(0, 0.5, 1, 3, Inf),
                     labels = c("A", "B", "C", "D"))
toy_dat <- cbind(toy_dat, as.data.frame(x_mat[, - 1L, drop = FALSE]))

## model-fitting from given model formula
fit2 <- cox_cure(~ x2 + x3 + x4 + group, ~ x2 + x3 + group,
                 time = time, event = status, data = toy_dat,
                 subset = group != "D", bootstrap = 30)
summary(fit2)

## get BIC's
BIC(fit1)
BIC(fit2)
BIC(fit1, fit2)


### Cox cure rate model with uncertain event status ==================
## simulate sample data
sim_dat <- simData4cure(nSubject = 200, max_censor = 10, lambda_censor = 0.1,
                        survMat = x_mat, cureMat = x_mat, b0 = 1)
table(sim_dat$case)
table(sim_dat$obs_event, useNA = "ifany")

## use formula
fit3 <- cox_cure(~ x1 + x2 + x3, ~ z1 + z2 + z3,
                 time = obs_time, event = obs_event, data = sim_dat)
summary(fit3)

## use design matrix
fit4 <- cox_cure.fit(x_mat, x_mat, time = sim_dat$obs_time,
                     event = sim_dat$obs_event)
summary(fit4)

## get BIC's
BIC(fit3, fit4)
