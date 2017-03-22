## estimate se by bootstrapping methods
bootSe <- function(obj, numBoot = 200, control = list(), ...)
{
    fm <- obj@formula
    cal <- obj@call
    censorRate0 <- obj@start$censorRate0
    control <- do.call(control_bootSe,
                       c(control, list(censorRate0_ = censorRate0)))
    start0 <- list(beta = obj@start$beta, censorRate = control$startGrid)
    cal$start <- quote(start0)
    cal$control <- obj@control
    cal$control$noSE_ <- TRUE
    cal$data <- quote(bootDat)
    dat <- obj@data
    idName <- as.character(fm[[2L]][[2L]])
    timeName <- as.character(fm[[2L]][[3L]])
    dat <- dat[order(dat[, idName], dat[, timeName]), ]
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
    if (control$estOnly)
        return(as.numeric(estMat))
    obj@estimates$beta[, "se(coef)"] <- apply(estMat, 1L, sd)
    tmp <- obj@estimates$beta[, "z"] <- obj@estimates$beta[, "coef"] /
        obj@estimates$beta[, "se(coef)"]
    obj@estimates$beta[, "Pr(>|z|)"] <- 2 * stats::pnorm(- abs(tmp))
    obj
}


### internal functions
control_bootSe <- function(startGrid, fixStart = FALSE, estOnly = FALSE,
                           ..., censorRate0_)
{
    startGrid <- if (fixStart) {
                     censorRate0_
                 } else if (missing(startGrid)){
                     seq.int(max(0, censorRate0_ - 0.2),
                             min(1, censorRate0_ + 0.2), 0.05)
                 }
    list(startGrid = startGrid, estOnly = estOnly)
}
