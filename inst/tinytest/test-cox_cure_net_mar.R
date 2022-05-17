library(intsurv)

## simulate a toy data
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
                    b0 = 0.5, p1 = 0.95, p2 = 0.95, p3 = 0.95)

## model-fitting from given design matrices
## MCAR
fit11 <- cox_cure_net.fit(
    x_mat,
    x_mat,
    dat$obs_time,
    dat$obs_event,
    surv_control = list(nlambda = 5, alpha = 0.8),
    cure_control = list(nlambda = 5, alpha = 0.8)
)

## MAR


## model-fitting from given model formula
fm <- paste(paste0("x", seq_len(p)), collapse = " + ")
surv_fm <- as.formula(sprintf("~ %s", fm))
cure_fm <- surv_fm
fit21 <- cox_cure_net(surv_fm, cure_fm, data = dat,
                     time = obs_time, event = obs_event,
                     surv_nlambda = 5, cure_nlambda = 5,
                     surv_alpha = 0.5, cure_alpha = 0.5)

## summary of BIC's
BIC(fit11)
BIC(fit21)

## list of coefficient estimates based on BIC
coef(fit1)
coef(fit2)
