library(intsurv)
## generate simulated survival data with uncertain records
set.seed(123)
simuDat <- simData4iCoxph(nSubject = 200)
## fit the integertive Cox model
fit <- iCoxph(Survi(ID, time, event) ~ x1 + x2 + x3 + x4,
              data = simuDat, start = iCoxph.start(methods = "nearest"),
              control = iCoxph.control(tol_beta = 1e-5))
## estimated covariate coefficients
coef(fit)
## get SE estimates by bootstrap
fit <- bootSe(fit, B = 30)
## summary of the fitted model
summary(fit)
