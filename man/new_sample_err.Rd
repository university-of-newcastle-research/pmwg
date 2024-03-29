% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/messaging.R
\name{new_sample_err}
\alias{new_sample_err}
\title{Error handler forany error in new_sample function call(s)}
\usage{
new_sample_err(pmwgs, envir, err_cond)
}
\arguments{
\item{pmwgs}{The pmwgs object for the current run.}

\item{envir}{The environment of the function at this point in time.}

\item{err_cond}{The original error condition that prompted this.}
}
\description{
If an error was detected when generating new samples. Save the state
of the samples and particles at that moment to help with debugging.
}
\keyword{internal}
