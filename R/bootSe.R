## estimate se by bootstrapping methods
bootSe <- function(obj, numBoot = 50, fixStart = FALSE, ...) {
    cal <- obj@call
    censorRate0 <- obj@start$censorRate0
    start0 <- if (fixStart) {
                  ## get the starting value used
                  list(beta = obj@start$beta, censorRate = censorRate0)
              } else {
                  seq_cen <- seq.int(max(0, censorRate0 - 0.2),
                                     min(1, censorRate0 + 0.2), 0.05)
                  list(beta = obj@start$beta,
                       censorRate = seq_cen)
              }
    cal$start <- quote(start0)
    cal$control <- obj@control
    cal$control$noSE_ <- TRUE
    cal$data <- quote(bootDat)
    dat <- obj@data
    dat <- dat[with(dat, order(ID, obsTime)), ]
    idTab <- table(dat$ID)
    uid <- unique(dat$ID)
    nSub <- length(uid)
    estMat <- replicate(numBoot, {
        sID <- sort(sample(uid, size = nSub, replace = TRUE))
        tmpDat <- data.frame(ID = sID)
        repNum <- idTab[match(as.character(sID), names(idTab))]
        bootDat <- merge(tmpDat, dat, by = "ID")
        bootDat$ID <- rep(seq_along(uid), repNum)
        res <- eval(cal)
        as.numeric(res@estimates$beta[, "coef"])
    })
    seEst <- apply(estMat, 1L, sd)
    est <- rowMeans(estMat)
    list(se = seEst, beta = est)
}
