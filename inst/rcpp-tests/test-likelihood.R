library(tinytest)
Rcpp::sourceCpp("test-likelihood.cpp")

## simulate right-censored data with a cure fraction
set.seed(123)
n_obs <- 1e3
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
rord <- sample(n_obs)
train_surv_offset <- runif(n_obs)
train_cure_offset <- runif(n_obs)
train_mar_offset <- runif(n_obs)

## test cox cure rate model: zero tail completion
res1 <- rt_ll_CoxphCure(
    train_surv_x = x_mat,
    train_cure_x = x_mat[, - 1],
    train_time = obs_time,
    train_event = event,
    test_surv_x = x_mat[rord, ],
    test_cure_x = x_mat[rord, - 1],
    test_time = obs_time[rord],
    test_event = event[rord],
    train_surv_offset = train_surv_offset,
    train_cure_offset = train_cure_offset,
    test_surv_offset = train_surv_offset[rord],
    test_cure_offset = train_cure_offset[rord],
    surv_standardize = FALSE,
    cure_standardize = TRUE,
    tail_completion = 1
)

expect_equal(res1$train_negLogL, res1$test_negLogL, tolerance = 1e-5)
expect_equal(as.numeric(res1$new_time), obs_time[ord])
expect_equal(as.numeric(res1$new_event), event[ord])
expect_false(res1$new_surv_standaridze)
expect_true(res1$new_cure_standaridze)
expect_equal(c(res1$ord), order(obs_time[rord], - event[rord]) - 1L)

rev_ord <- order(c(res1$ord))
expect_equal(res1$new_surv_offset[rev_ord], train_surv_offset[rord])
expect_equal(res1$new_cure_offset[rev_ord], train_cure_offset[rord])

## test cox cure rate model: exp tail completion
res2 <- rt_ll_CoxphCure(
    train_surv_x = x_mat,
    train_cure_x = x_mat[, - 1],
    train_time = obs_time,
    train_event = event,
    test_surv_x = x_mat[rord, ],
    test_cure_x = x_mat[rord, - 1],
    test_time = obs_time[rord],
    test_event = event[rord],
    train_surv_offset = train_surv_offset,
    train_cure_offset = train_cure_offset,
    test_surv_offset = train_surv_offset[rord],
    test_cure_offset = train_cure_offset[rord],
    surv_standardize = FALSE,
    cure_standardize = TRUE,
    tail_completion = 2
)

expect_equal(res2$train_negLogL, res2$test_negLogL, tolerance = 1e-5)
expect_equal(as.numeric(res2$new_time), obs_time[ord])
expect_equal(as.numeric(res2$new_event), event[ord])
expect_false(res2$new_surv_standaridze)
expect_true(res2$new_cure_standaridze)
expect_equal(c(res2$ord), order(obs_time[rord], - event[rord]) - 1L)

rev_ord <- order(c(res2$ord))
expect_equal(res2$new_surv_offset[rev_ord], train_surv_offset[rord])
expect_equal(res2$new_cure_offset[rev_ord], train_cure_offset[rord])


## test cox cure rate mar: zero tail
event2 <- event
event2[sample(n_obs, size = 100)] <- 0.5

res3 <- rt_ll_CoxphCureMar(
    train_surv_x = x_mat,
    train_cure_x = x_mat[, - 1],
    train_mar_x = x_mat[, - 2],
    train_time = obs_time,
    train_event = event2,
    test_surv_x = x_mat[rord, ],
    test_cure_x = x_mat[rord, - 1],
    test_mar_x = x_mat[rord, - 2],
    test_time = obs_time[rord],
    test_event = event2[rord],
    train_surv_offset = train_surv_offset,
    train_cure_offset = train_cure_offset,
    train_mar_offset = train_mar_offset,
    test_surv_offset = train_surv_offset[rord],
    test_cure_offset = train_cure_offset[rord],
    test_mar_offset = train_mar_offset[rord],
    surv_standardize = FALSE,
    cure_standardize = TRUE,
    tail_completion = 1
)

expect_equal(res3$train_negLogL, res3$test_negLogL, tolerance = 1e-5)

ord <- order(obs_time, - event2)
expect_equal(as.numeric(res3$new_time), obs_time[ord])
expect_equal(as.numeric(res3$new_event), event2[ord])
expect_false(res3$new_surv_standaridze)
expect_true(res3$new_cure_standaridze)
expect_equal(c(res3$ord), order(obs_time[rord], - event2[rord]) - 1L)

rev_ord <- order(c(res3$ord))
expect_equal(res3$new_surv_offset[rev_ord], train_surv_offset[rord])
expect_equal(res3$new_cure_offset[rev_ord], train_cure_offset[rord])
expect_equal(res3$new_mar_offset[rev_ord], train_mar_offset[rord])

## test cox cure rate mar: exp tail
res4 <- rt_ll_CoxphCureMar(
    train_surv_x = x_mat,
    train_cure_x = x_mat[, - 1],
    train_mar_x = x_mat[, - 2],
    train_time = obs_time,
    train_event = event2,
    test_surv_x = x_mat[rord, ],
    test_cure_x = x_mat[rord, - 1],
    test_mar_x = x_mat[rord, - 2],
    test_time = obs_time[rord],
    test_event = event2[rord],
    train_surv_offset = train_surv_offset,
    train_cure_offset = train_cure_offset,
    train_mar_offset = train_mar_offset,
    test_surv_offset = train_surv_offset[rord],
    test_cure_offset = train_cure_offset[rord],
    test_mar_offset = train_mar_offset[rord],
    surv_standardize = FALSE,
    cure_standardize = TRUE,
    tail_completion = 2
)

expect_equal(res4$train_negLogL, res4$test_negLogL, tolerance = 1e-5)

ord <- order(obs_time, - event2)
expect_equal(as.numeric(res4$new_time), obs_time[ord])
expect_equal(as.numeric(res4$new_event), event2[ord])
expect_false(res4$new_surv_standaridze)
expect_true(res4$new_cure_standaridze)
expect_equal(c(res4$ord), order(obs_time[rord], - event2[rord]) - 1L)

rev_ord <- order(c(res4$ord))
expect_equal(res4$new_surv_offset[rev_ord], train_surv_offset[rord])
expect_equal(res4$new_cure_offset[rev_ord], train_cure_offset[rord])
expect_equal(res4$new_mar_offset[rev_ord], train_mar_offset[rord])
