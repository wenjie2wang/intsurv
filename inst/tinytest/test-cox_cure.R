library(intsurv)

## create a new environment for evaluation
## env_a <- new.env(parent = .GlobalEnv)
env_a <- new.env()
env_b <- new.env(parent = env_a)

with(env_b, {
    ## simulate a dataset
    n_obs <- 2e2
    p <- 5
    x_mat <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
    colnames(x_mat) <- paste0("X", seq_len(p))
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
    sim_dat <- data.frame(time = obs_time, status = event)
    sim_dat <- cbind(sim_dat, x_mat[, - 1L])
    sim_dat$group <- cut(abs(x_mat[, 1L]), breaks = c(0, 0.5, 1, 3, Inf),
                         labels = c("A", "B", "C", "D"))
    ## test the formula interface with data
    fit <- cox_cure(~ X2 + X3 + X4 + group,
                    ~ X2 + X3 + group,
                    time = time,
                    event = status,
                    data = sim_dat,
                    subset = group != "D",
                    bootstrap = 10)
    fit
    summary(fit)
})

with(env_b, {
    ## test the formula interface without data
    fit <- with(
        sim_dat,
        cox_cure(~ X2 + X3 + X4 + group,
                 ~ X2 + X3 + group,
                 time = time,
                 event = status,
                 subset = group != "D",
                 bootstrap = 10)
    )
    summary(fit)
    ## test the matrix interface
    (fit <- cox_cure.fit(x_mat, x_mat, obs_time, event, bootstrap = 10))
    summary(fit)
})

## expect errors from prep_cure_model
with(env_b, {
    ## surv formula missing
    expect_error(
        cox_cure(cure_formula = ~ X2 + X3),
        pattern = "formula"
    )
    ## cure formula missing
    expect_error(
        cox_cure(~ X2 + X3),
        pattern = "formula"
    )
    ## time missing
    expect_error(
        cox_cure(~ X2 + X3, ~ X2 + X3, data = sim_dat, event = status),
        pattern = "time"
    )
    ## event indicators missing
    expect_error(
        cox_cure(~ X2 + X3, ~ X2 + X3, data = sim_dat, time = time),
        pattern = "event"
    )
})

## check backward compatibility
## simulate a dataset
n_obs <- 2e2
p <- 5
x_mat <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
colnames(x_mat) <- paste0("X", seq_len(p))
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
sim_dat <- data.frame(time = obs_time, status = event)
sim_dat <- cbind(sim_dat, x_mat)
## test the formula interface with data
fm <- paste(paste0("X", seq_len(p)), collapse = " + ")
surv_fm <- as.formula(sprintf("~ %s", fm))
cure_fm <- surv_fm

## em_max_iter, em_rel_tol, tail_completion
expect_warning(
    old_fit <- cox_cure.fit(x_mat, x_mat, obs_time, event,
                            em_max_iter = 3,
                            em_rel_tol = 1e-3,
                            tail_completion = 2),
    pattern = "cox_cure.fit()"
)
expect_equal(old_fit$control$maxit, 3)
expect_equal(old_fit$control$epsilon, 1e-3)
expect_equal(old_fit$control$tail_completion, 2)

## firth, ealry_stop
expect_warning(
    cox_cure.fit(x_mat, x_mat, obs_time, event, firth = TRUE),
    pattern = "firth"
)
expect_warning(
    cox_cure.fit(x_mat, x_mat, obs_time, event, early_stop = TRUE),
    pattern = "early_stop"
)

## surv_start, cure_start
start0 <- rep(0.1, ncol(x_mat))
expect_warning(
    old_fit <- cox_cure.fit(x_mat, x_mat, obs_time, event,
                            surv_start = start0,
                            cure_start = start0),
    pattern = "surv_start"
)
expect_equal(old_fit$surv_mstep$start, start0)
expect_equal(old_fit$cure_mstep$start, start0)

## surv_offset, cure_offset
set.seed(1)
offset0 <- runif(nrow(x_mat))
expect_warning(
    old_fit <- cox_cure.fit(x_mat, x_mat, obs_time, event,
                            surv_offset = offset0,
                            cure_offset = offset0),
    pattern = "surv_offset"
)
expect_equal(old_fit$surv_mstep$offset, offset0)
expect_equal(old_fit$cure_mstep$offset, offset0)

## surv_max_iter, cure_max_iter
max_iter0 <- 3L
expect_warning(
    old_fit <- cox_cure.fit(x_mat, x_mat, obs_time, event,
                            surv_max_iter = max_iter0,
                            cure_max_iter = max_iter0),
    pattern = "surv_max_iter"
)
expect_equal(old_fit$surv_mstep$maxit, max_iter0)
expect_equal(old_fit$cure_mstep$maxit, max_iter0)
