##' For right-censored data with uncertain event status, fit the Cox cure model
##' proposed by Wang et al. (2020).
##'
##' @param event   \code{NA}'s are
##'     allowed and represent uncertain event status.
##'
##' @return A \code{cox_cure_mar} object for the fitted Cox cure rate model with
##'     uncertain event indicators.
##'
##' @references
##' Wang, W., Luo, C., Aseltine, R. H., Wang, F., Yan, J., & Chen, K. (2020).
##' Suicide Risk Modeling with Uncertain Diagnostic Records. \emph{arXiv
##' preprint arXiv:2009.02597}.
##'
##' @export
cox_cure_mar <- function(surv_formula,
                         cure_formula,
                         mar_forumla,
                         time,
                         event,
                         data,
                         subset,
                         contrasts = NULL,
                         bootstrap = 0,
                         control = cox_cure.control(),
                         surv_control = intsurv.control(),
                         cure_control = intsurv.control(),
                         mar_control = intsurv.control(),
                         ...)
{

}
