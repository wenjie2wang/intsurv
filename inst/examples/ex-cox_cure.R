library(intsurv)

### ordinary Cox cure rate model =====================================
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
                     labels = LETTERS[1:4])
toy_dat <- cbind(toy_dat, as.data.frame(x_mat[, - 1L, drop = FALSE]))

## model-fitting from given model formula
fit2 <- cox_cure(
    ~ x3 + x4 + group,
    ~ group + x3 + offset(x2),
    time = time,
    event = status,
    data = toy_dat,
    subset = group != "D",
    bootstrap = 30
)
summary(fit2)

## get BIC's
BIC(fit1)
BIC(fit2)
BIC(fit1, fit2)


### regularized Cox cure rate model ==================================
## simulate a toy right-censored data with a cure fraction
set.seed(123)
n_obs <- 100
p <- 10
x_mat <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
colnames(x_mat) <- paste0("x", seq_len(p))
surv_beta <- c(rep(0, p - 5), rep(1, 5))
cure_beta <- c(rep(1, 2), rep(0, p - 2))
dat <- simData4cure(nSubject = n_obs, lambda_censor = 0.01,
                    max_censor = 10, survMat = x_mat,
                    survCoef = surv_beta, cureCoef = cure_beta,
                    b0 = 0.5, p1 = 1, p2 = 1, p3 = 1)

## model-fitting from given design matrices
fit1 <- cox_cure_net.fit(x_mat, x_mat, dat$obs_time, dat$obs_event,
                         surv_control = list(nlambda = 10, alpha = 0.8),
                         cure_control = list(nlambda = 10, alpha = 0.8))

## model-fitting from given model formula
fm <- paste(paste0("x", seq_len(p)), collapse = " + ")
surv_fm <- as.formula(sprintf("~ %s", fm))
cure_fm <- surv_fm
fit2 <- cox_cure_net(surv_fm, cure_fm, data = dat,
                     time = obs_time, event = obs_event)

## summary of BIC's
BIC(fit1)
BIC(fit2)
BIC(fit1)[which.min(BIC(fit1)[, "BIC"]), ]
BIC(fit2)[which.min(BIC(fit2)[, "BIC"]), ]

## list of coefficient estimates based on BIC
coef(fit1)
coef(fit2)
