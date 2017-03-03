## estimate se by bootstrapping methods
bootSe <- function(obj, numBoot = 50, fixStart = FALSE, ...) {
    cal <- obj@call
    if (fixStart) {
        ## get the starting value used
        censorRate0 <- obj@start$censorRate0
        start0 <- list(beta = obj@start$beta, censorRate = censorRate0)
        cal$start <- quote(start0)
    }
    cal$data <- quote(bootDat)
    dat <- obj@data
    uid <- unique(dat$ID)
    nSub <- length(uid)
    estMat <- replicate(numBoot, {
        sID <- sample(uid, size = nSub, replace = TRUE)
        bootDat <- base::subset(dat, ID %in% sID)
        res <- eval(cal)
        as.numeric(res@estimates$beta[, "coef"])
    })
    apply(estMat, 1L, sd)
}
