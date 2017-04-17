##' Formula Response for Survival Data With Uncertainty
##'
##' \code{survi} returns an S4 class that represents formula response for
##' survival data with uncertain records due to imperfect data integration.  The
##' last letter 'i' in 'survi' represents 'integration'.
##'
##' @param ID Identificator of each subject.
##' @param time Time of reccurence event or censoring.
##' @param event The status indicator, 0 = censored, 1 = event.
##' @aliases survi
##' @export
survi <- function(ID, time, event)
{
    ## inpDat <- data.frame(id, time, event)
    ## TODO: function to check the date structure of uncertained records
    ## dat <- checkSurvi(inpDat)
    ## outDat <- with(dat, as.matrix(cbind(ID, time, event)))
    mat <- as.matrix(cbind(ID, time, event))
    methods::new("survi", mat,
                 ID = as.character(ID),
                 time = as.numeric(time),
                 event = as.integer(event))
}
