### prepare design matrix and response from formula inputs
prep_model <- function(formula, time, event, x, data,
                       subset, na.action, contrasts = NULL)
{
    this_call <- match.call(expand.dots = FALSE)
    ## if formula missing, use given x
    if (missing(formula)) {
        if (missing(x)) {
            stop("Neither the design matrix 'x' nor the 'formula' is missing.",
                 call. = FALSE)
        }
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
    mm <- stats::model.matrix.default(formula, mf,
                                      contrasts.arg = contrasts)
    ## output: na.action
    na.action <- if (is.null(attr(mf, "na.action"))) {
                     options("na.action")[[1L]]
                 } else {
                     paste0("na.", class(attr(mf, "na.action")))
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
