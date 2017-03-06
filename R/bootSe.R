## estimate se by bootstrapping methods
bootSe <- function(obj, numBoot = 50, fixStart = FALSE, ...) {
    cal <- obj@call
    if (fixStart) {
        ## get the starting value used
        censorRate0 <- obj@start$censorRate0
        start0 <- list(beta = obj@start$beta, censorRate = censorRate0)
        cal$start <- quote(start0)
    } else {
        censorRate0 <- obj@start$censorRate0
        seq_cen <- seq.int(max(0, censorRate0 - 0.2),
                           min(1, censorRate0 + 0.2), 0.02)
        start0 <- list(beta = obj@start$beta,
                       censorRate = seq_cen)
        cal$start <- quote(start0)
    }
    cal$data <- quote(bootDat)
    dat <- obj@data
    uid <- unique(dat$ID)
    nSub <- length(uid)
    estMat <- replicate(numBoot, {
        sID <- sample(uid, size = nSub, replace = TRUE)
        tmpDat <- data.frame(ID = sID)
        bootDat <- merge(tmpDat, dat, by = "ID")
        res <- eval(cal)
        as.numeric(res@estimates$beta[, "coef"])
    })
    seEst <- apply(estMat, 1L, sd)
    est <- rowMeans(estMat)
    list(se = seEst, beta = est)
}
