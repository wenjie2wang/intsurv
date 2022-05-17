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
## MCAR
fit31 <- cox_cure(
    ~ x1 + x2 + x3,
    ~ z1 + z2 + z3,
    time = obs_time,
    event = obs_event,
    data = sim_dat
)
summary(fit31)

## MAR
fit32 <- cox_cure_mar(
    ~ x1 + x2 + x3,
    ~ z1 + z2 + z3,
    ~ x1 + x2,
    time = obs_time,
    event = obs_event,
    data = sim_dat
)
summary(fit32)

## use design matrix
## MCAR
fit41 <- cox_cure.fit(
    x_mat,
    x_mat,
    time = sim_dat$obs_time,
    event = sim_dat$obs_event
)
summary(fit41)

## MAR
fit42 <- cox_cure_mar.fit(
    x_mat,
    x_mat,
    x_mat,
    time = sim_dat$obs_time,
    event = sim_dat$obs_event
)
summary(fit42)

## get BIC's
BIC(fit31, fit41)
BIC(fit32, fit42)
