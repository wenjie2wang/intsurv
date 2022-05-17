library(intsurv)

### Cox cure rate model for uncertain event status ===================
set.seed(123)
n_obs <- 200
p <- 5
x_mat <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
colnames(x_mat) <- paste0("x", seq_len(p))

## simulate sample data
sim_dat <- simData4cure(nSubject = 200, max_censor = 10, lambda_censor = 0.1,
                        survMat = x_mat, cureMat = x_mat, b0 = 1)
table(sim_dat$case)
table(sim_dat$obs_event, useNA = "ifany")

## use formula
fit3 <- cox_cure_mar(
    ~ x1 + x2 + x3,
    ~ z1 + z2 + z3,
    ~ x1 + x2 + x3,
    time = obs_time,
    event = obs_event,
    data = sim_dat
)
summary(fit3)

## use design matrix
fit4 <- cox_cure_mar.fit(x_mat, x_mat, x_mat,
                         time = sim_dat$obs_time,
                         event = sim_dat$obs_event)
summary(fit4)

## get BIC's
BIC(fit3, fit4)

### regularized Cox cure model for uncertain event status ============
## simulate a toy data
n_obs <- 100
p <- 10
x_mat <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
colnames(x_mat) <- paste0("x", seq_len(p))
surv_beta <- c(rep(0, p - 5), rep(1, 5))
cure_beta <- c(rep(1, 2), rep(0, p - 2))
dat <- simData4cure(nSubject = n_obs, lambda_censor = 0.01,
                    max_censor = 10, survMat = x_mat,
                    survCoef = surv_beta, cureCoef = cure_beta,
                    b0 = 0.5, p1 = 0.95, p2 = 0.95, p3 = 0.95)

## model-fitting from given design matrices
fit1 <- cox_cure_net_mar.fit(
    x_mat, x_mat, x_mat,
    dat$obs_time, dat$obs_event,
    surv_control = list(nlambda = 5, alpha = 0.8)
)

## model-fitting from given model formula
fm <- paste(paste0("x", seq_len(p)), collapse = " + ")
surv_fm <- as.formula(sprintf("~ %s", fm))
cure_fm <- surv_fm
fit2 <- cox_cure_net_mar(
    surv_fm,
    cure_fm,
    cure_fm,
    data = dat,
    time = obs_time,
    event = obs_event,
    cure_control = list(nlambda = 5, alpha = 1),
    mar_control = list(nlambda = 10, alpha = 0.9)
)

## summary of BIC's
BIC(fit1)
BIC(fit2)

## list of coefficient estimates based on BIC
coef(fit1)
coef(fit2)
