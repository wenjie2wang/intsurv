##' Formula Response for Survival Data With Uncertainty
##'
##' \code{survi} returns an S4 class that represents formula response for
##' survival data with uncertain records due to imperfect data integration.  The
##' last letter 'i' in 'survi' represents 'integration'.
##'
##' @param id Identificator of each subject.
##' @param time Time of reccurence event or censoring.
##' @param event The status indicator, 0 = censored, 1 = event.
##' @aliases survi
##' @export
survi <- function(id, time, event)
{
    ## inpDat <- data.frame(id, time, event)
    ## TODO: function to check the date structure of uncertained records
    ## dat <- checkSurvi(inpDat)
    ## outDat <- with(dat, as.matrix(cbind(id, time, event)))
    mat <- as.matrix(cbind(id, time, event))
    methods::new("survi", mat,
                 id = as.character(id),
                 time = as.numeric(time),
                 event = as.integer(event))
}
