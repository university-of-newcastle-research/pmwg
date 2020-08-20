#' pmwg: Particle Metropolis Within Gibbs.
#'
#' The pmwg package provides a general purpose implementation of the
#' sampling techniques outlined in
#' \href{https://doi.org/10.1016/j.jmp.2020.102368}{Gunawan et al. (2020)}.
#' The user of this package is required to provide their own log likelihood
#' function, but given this the functions provided can estimate model
#' parameters, the full covariance matrix and subject random effects in a
#' hierarchical Bayesian way.
#'
#' @section Documentation:
#' The documentation found at \url{https://newcastlecl.github.io/samplerDoc/}
#' contains background information and motivation for the approach used in
#' this package and several detailed examples of the package in action. It also
#' includes a list of common problems and associated troubleshooting steps.
#'
#' @section User input:
#' The user is expected to provide a data source in a format that is compatible
#' with R data.frame methods. This data must have at least one column named
#' `subject` that has a unique identifier for each participants data.
#'
#' Additionally the user should provide a function that when given a set of
#' parameter estimates and the data for a single subject return the log of the
#' likelihood of that data given the parameter estimates.
#'
#' The final piece of required information is a list of the names of each
#' parameter that should be estimated. There is also the capability to provide
#' priors on the model parameters, start points for the model parameters and
#' covariance matrix as well as options to fine tune the sampling process
#'
#' @references
#' Gunawan, D., Hawkins, G. E., Tran, M. N., Kohn, R., & Brown, S. D. (2020).
#' New estimation approaches for the hierarchical Linear Ballistic Accumulator
#' model. \emph{Journal of Mathematical Psychology, 96}, 102368.
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
