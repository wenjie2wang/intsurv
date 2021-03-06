% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootSe.R
\name{bootSe}
\alias{bootSe}
\title{Standard Error Estimates through Bootstrap}
\usage{
bootSe(
  object,
  B = 50,
  se = c("inter-quartile", "mad", "sd"),
  return_beta = FALSE,
  ...
)
}
\arguments{
\item{object}{\code{\link{iCoxph-class}} object.}

\item{B}{A positive integer specifying number of bootstrap samples used for
SE estimates.  A large number, such as 200, is often needed for a more
reliable estimation in practice.  If \code{B = 1} is specified, the
function will return the covariate coefficient estimates instead of a
\code{iCoxph-class} object.}

\item{se}{A character value specifying the way computing SE from bootstrap
samples. The default method is based on median absolute deviation and
the second method is based on inter-quartile, both of which are based on
normality of the bootstrap estimates and provides robust estimates for
SE. The third method estimates SE by the standard deviation of the
bootstrap estimates.}

\item{return_beta}{A logical value. If \code{TRUE}, the function returns the
covariate coefficient estimates from the given number of bootstrap
samples, which allows users to split the most computationally intensive
step into small pieces that can be computed in a parallel manner.  The
default value is \code{FALSE}.}

\item{...}{Other arguments for future usage.  A warning will be thrown if
any invalid argument is specified.}
}
\value{
\code{\link{iCoxph-class}} object or a numeric matrix that contains
    the covariate coefficient estimates from the given number of bootstrap
    samples in rows.
}
\description{
For \code{\link{iCoxph-class}} object, add (or update) standard error (SE)
estimates through bootstrap methods, or compute the coefficient estimates
from the given number of bootstrap samples.
}
\details{
Three different methods are available for computing SE from bootstrap
samples through argument \code{se}.  Given the fact that the bootstrap
method is computationally intensive, the function returns the coefficient
estimates in a matrix from the given number of bootstrap samples when
\code{return_beta = TRUE)} is specified, which can be used in parallel
computing or high performance computing (HPC) cluster.  The SE estimates can
be further computed based on estimates from bootstrap samples by users on
their own.  The \code{return_beta = TRUE} is implied, when \code{B = 1} is
specified.
}
\examples{
## See examples of function 'iCoxph'.
}
\seealso{
\code{\link{iCoxph}} for fitting integrative Cox model.
}
