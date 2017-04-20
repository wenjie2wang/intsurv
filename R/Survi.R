################################################################################
##
##   R package intsurv by Wenjie Wang, Kun Chen and Jun Yan
##   Copyright (C) 2017
##
##   This file is part of the R package intsurv.
##
##   The R package intsurv is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package intsurv is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################


## collation after class.R
##' @include class.R
NULL


##' Formula Response for Survival Data With Uncertainty
##'
##' \code{Survi} returns an S4 class that represents formula response for
##' survival data with uncertain records due to imperfect data integration.  The
##' last letter 'i' in \code{Survi} represents 'integration'.
##'
##' @param ID Identificator of each subject.
##' @param time Time of reccurence event or censoring.
##' @param event The status indicator, 0 = censored, 1 = event.
##' @aliases Survi
##' @export
Survi <- function(ID, time, event)
{
    ## inpDat <- data.frame(id, time, event)
    ## TODO: function to check the date structure of uncertained records
    ## dat <- checkSurvi(inpDat)
    ## outDat <- with(dat, as.matrix(cbind(ID, time, event)))
    mat <- as.matrix(cbind(ID, time, event))
    methods::new("Survi", mat,
                 ID = as.character(ID),
                 time = as.numeric(time),
                 event = as.integer(event))
}
