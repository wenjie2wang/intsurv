library(intsurv)

## 1. simulate right-censored data with a cure fraction
set.seed(123)
n_obs <- 2e2
p <- 5
x_mat <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
colnames(x_mat) <- paste0("x", seq_len(p))
cure_beta <- rep(0.5, p)
b0 <- - 1
expit <- stats::binomial()$linkinv
ncure_prob <- expit(as.numeric(b0 + x_mat %*% cure_beta))
is_cure <- 1 - stats::rbinom(n_obs, size = 1, prob = ncure_prob)
surv_beta <- rep(0.5, p)
risk_score <- as.numeric(x_mat %*% surv_beta)
event_time <- rexp(n_obs, exp(as.numeric(x_mat %*% surv_beta)))
censor_time <- 10
event <- ifelse(event_time < censor_time & ! is_cure, 1, 0)
obs_time <- ifelse(event > 0, event_time, censor_time)

## model-fitting from given design matrices
cureFit1 <- cox_cure.fit(x_mat, x_mat, obs_time, event, bootstrap = 30)

## 2. create a toy dataset
toy_dat <- data.frame(time = obs_time, status = event)
toy_dat$group <- cut(abs(x_mat[, 1L]), breaks = c(0, 0.5, 1, 3, Inf),
                     labels = c("A", "B", "C", "D"))
toy_dat <- cbind(toy_dat, as.data.frame(x_mat[, - 1L, drop = FALSE]))

## model-fitting from given model formula
cureFit2 <- cox_cure(~ x2 + x3 + x4 + group, ~ x2 + x3 + group,
                     time = time, event = status, data = toy_dat,
                     subset = group != "D", bootstrap = 30)
