### prepare design matrix and response from formula inputs
prep_model <- function(formula, time, event, data,
                       subset, na.action, contrasts = NULL)
{
    this_call <- match.call(expand.dots = FALSE)
    if (missing(formula)) {
        stop("The 'formula' cannot be missing.",
             call. = FALSE)
    }
    ## if formula is a design matrix
    if (idx1 <- is.matrix(formula)) {
        x_names <- colnames(formula)
        formula <- as.formula(substitute(~ x))
        environment(formula) <- parent.frame()
        this_call$formula <- formula
    }
    if (missing(data)) {
        this_call$data <- environment(formula)
    }
    matched_call <- match(c("time", "event", "formula", "data",
                            "subset", "na.action"),
                          names(this_call), nomatch = 0L)
    this_call <- this_call[c(1L, matched_call)]

    ## drop unused levels in factors
    this_call$drop.unused.levels <- TRUE
    this_call[[1L]] <- quote(stats::model.frame.default)
    mf <- eval(this_call, parent.frame())
    mt <- attr(mf, "terms")
    ## suppress warnings on not used contrasts
    suppressWarnings({
        mm <- stats::model.matrix.default(formula, mf,
                                          contrasts.arg = contrasts)
    })
    ## output: na.action
    na.action <- if (is.null(attr(mf, "na.action"))) {
                     options("na.action")[[1L]]
                 } else {
                     paste0("na.", class(attr(mf, "na.action")))
                 }
    ## use original names
    if (idx1) {
        colnames(mm)[- 1L] <- x_names
    }
    ## output: contrasts
    contrasts <- attr(mm, "contrasts")
    ## return a list
    list(time = mf[["(time)"]],
         event = mf[["(event)"]],
         xMat = mm,
         na.action = na.action,
         contrasts = contrasts,
         xlevels = .getXlevels(mt, mf))
}
