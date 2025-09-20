library(intsurv)

### cox_cure for uncertain event status
n_obs <- 2e2
p <- 5
x_mat <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
colnames(x_mat) <- paste0("X", seq_len(p))
## simulate sample data
sim_dat <- simData4cure(nSubject = 200, max_censor = 10, lambda_censor = 0.1,
                        survMat = x_mat, cureMat = x_mat, b0 = 1)
table(sim_dat$case)
table(sim_dat$obs_event, useNA = "ifany")

## use formula
fit31 <- cox_cure(
    ~ x1 + x2 + x3,
    ~ z1 + z2 + z3,
    time = obs_time,
    event = obs_event,
    data = sim_dat
)
summary(fit31)

## use design matrix
fit41 <- cox_cure.fit(
    x_mat,
    x_mat,
    time = sim_dat$obs_time,
    event = sim_dat$obs_event
)
summary(fit41)

## get BIC's
BIC(fit31, fit41)
