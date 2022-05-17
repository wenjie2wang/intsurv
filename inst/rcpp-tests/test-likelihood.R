library(tinytest)
Rcpp::sourceCpp("test-likelihood.cpp")

## simulate right-censored data with a cure fraction
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
ord <- order(obs_time, - event)

## test cox cure rate model
res1 <- rt_coxph_cure(
    train_surv_x = x_mat,
    train_cure_x = x_mat[, - 1],
    train_time = obs_time,
    train_event = event,
    test_surv_x = x_mat,
    test_cure_x = x_mat[, - 1],
    test_time = obs_time,
    test_event = event,
    surv_offset = runif(n_obs),
    cure_offset = runif(n_obs),
    surv_standardize = FALSE,
    cure_standardize = TRUE
)

expect_equal(res1$train_negLogL, res1$test_negLogL, tolerance = 1e-5)
expect_equivalent(res1$new_surv_x, x_mat[ord, ])
expect_equivalent(res1$new_cure_x, x_mat[ord, - 1])
expect_equal(as.numeric(res1$new_time), obs_time[ord])
expect_equal(as.numeric(res1$new_event), event[ord])

## test cox cure rate mar
res2 <- rt_coxph_cure_mar(
    train_surv_x = x_mat,
    train_cure_x = x_mat[, - 1],
    train_mar_x = x_mat[, - 2],
    train_time = obs_time,
    train_event = event,
    test_surv_x = x_mat,
    test_cure_x = x_mat[, - 1],
    test_mar_x = x_mat[, - 2],
    test_time = obs_time,
    test_event = event,
    surv_offset = runif(n_obs),
    cure_offset = runif(n_obs),
    mar_offset = runif(n_obs),
    surv_standardize = FALSE,
    cure_standardize = TRUE
)

expect_equal(res2$train_negLogL, res2$test_negLogL, tolerance = 1e-5)
expect_equivalent(res2$new_surv_x, x_mat[ord, ])
expect_equivalent(res2$new_cure_x, x_mat[ord, - 1])
expect_equal(as.numeric(res2$new_time), obs_time[ord])
expect_equal(as.numeric(res2$new_event), event[ord])
