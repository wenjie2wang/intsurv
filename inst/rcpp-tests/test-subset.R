library(tinytest)
Rcpp::sourceCpp("test-subset.cpp")

## simulate right-censored data with a cure fraction
set.seed(123)
n_obs <- 200
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
surv_offset <- runif(n_obs)
cure_offset <- runif(n_obs)
mar_offset <- runif(n_obs)
index <- sample(n_obs, size = 100)

## subset CoxphCure
res1 <- rt_subset_CoxphCure(
    index - 1,
    surv_x = x_mat,
    cure_x = x_mat[, - 1],
    time = obs_time,
    event = event,
    surv_offset = surv_offset,
    cure_offset = cure_offset,
    surv_standardize = TRUE,
    cure_standardize = FALSE
)

ord0 <- order(obs_time, - event)
sub_idx <- ord0[index]
expect_equal(obs_time[sub_idx], res1$new_time[res1$rev_ord + 1])
expect_equal(event[sub_idx], res1$new_event[res1$rev_ord + 1])
expect_equivalent(x_mat[sub_idx, ], res1$new_surv_x[res1$rev_ord + 1, ])
expect_equivalent(x_mat[sub_idx, - 1], res1$new_cure_x[res1$rev_ord + 1, ])
expect_equal(surv_offset[sub_idx], res1$new_surv_offset[res1$rev_ord + 1])
expect_equal(cure_offset[sub_idx], res1$new_cure_offset[res1$rev_ord + 1])
expect_true(res1$new_surv_standaridze)
expect_false(res1$new_cure_standaridze)

## subset CoxphCureMar
res2 <- rt_subset_CoxphCureMar(
    index - 1,
    surv_x = x_mat,
    cure_x = x_mat[, - 1],
    mar_x = x_mat[, - 2],
    time = obs_time,
    event = event,
    surv_offset = surv_offset,
    cure_offset = cure_offset,
    mar_offset = mar_offset,
    surv_standardize = TRUE,
    cure_standardize = FALSE,
    mar_standardize = TRUE
)

ord0 <- order(obs_time, - event)
sub_idx <- ord0[index]
expect_equal(obs_time[sub_idx], res2$new_time[res2$rev_ord + 1])
expect_equal(event[sub_idx], res2$new_event[res2$rev_ord + 1])
expect_equivalent(x_mat[sub_idx, ], res2$new_surv_x[res2$rev_ord + 1, ])
expect_equivalent(x_mat[sub_idx, - 1], res2$new_cure_x[res2$rev_ord + 1, ])
expect_equal(surv_offset[sub_idx], res2$new_surv_offset[res2$rev_ord + 1])
expect_equal(cure_offset[sub_idx], res2$new_cure_offset[res2$rev_ord + 1])
expect_equal(mar_offset[sub_idx], res2$new_mar_offset[res2$rev_ord + 1])
expect_true(res2$new_surv_standaridze)
expect_false(res2$new_cure_standaridze)
expect_true(res2$new_mar_standaridze)
