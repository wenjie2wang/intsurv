library(intsurv)

## simulate a dataset
n_obs <- 2e2
p <- 5
x_mat <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
colnames(x_mat) <- paste0("X", seq_len(p))
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
sim_dat <- data.frame(time = obs_time, status = event)
sim_dat <- cbind(sim_dat, x_mat[, - 1L])
sim_dat$group <- cut(abs(x_mat[, 1L]), breaks = c(0, 0.5, 1, 3, Inf),
                     labels = c("A", "B", "C", "D"))

## test the formula interface with data
fit <- cox_cure(~ X2 + X3 + X4 + group, ~ X2 + X3 + group,
                time = time, event = status,
                data = sim_dat, subset = group != "D",
                bootstrap = 10)
summary(fit)

## test the formula interface without data
fit <- with(
    sim_dat,
    cox_cure(~ X2 + X3 + X4 + group, ~ X2 + X3 + group,
             time = time, event = status,
             subset = group != "D",
             bootstrap = 10)
)
summary(fit)

## test the matrix interface
fit <- cox_cure.fit(x_mat, x_mat, obs_time, event, bootstrap = 10)
summary(fit)
