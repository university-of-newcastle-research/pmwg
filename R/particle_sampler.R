#' Initialise a PMwG sampler and return a pmwgs object
#'
#' Takes the initial parameter set and creates default values for a
#' range of PMwG necessary variables.
#'
#' @param data The data.frame containing empirical data to be modelled. Assumed
#'   to contain at least one column called subject whose elements are unique
#'   identifiers for each subject.
#' @param parameters The list of parameter names to be used in the model
#' @param burn A vector containing the number of particles and number of
#'   iterations respectively for the burnin phase.
#' @param adapt A vector containing the number of particles and maximum number
#'   of iterations for the adaptation phase.
#' @param sample A vector containing the number of particles and number of
#'   iterations respectively for the sampling phase.
#' @param llfunc A log likelihood function that given a list of parameter values
#'   and a data.frame containing subject data will return the log likelihood of
#'   \code{x} given \code{data}. See \code{\link{lba_loglike}} for an example.
#'
#' @return A pmwgs object that can be run, plotted and more
#' @examples
#' sampler <- pmwgs(forstmann, c("b1", "b2", "b3", "A", "v1", "v2", "t0"))
#' @export
pmwgs <- function(data, parameters, burn = c(1000, 500), adapt = c(1000, 5000),
                  sample = c(100, 1000), llfunc = lba_loglike) {
  n_pars <- length(parameters)
  # Tuning settings for the Gibbs steps
  hyper <- list(
    dof <- 2,  # hyperparameter on prior (Half-t degrees of freedom)
    scale <- 1 # hyperparameter on prior (Half-t scale)
  )
  n_subjects <- length(unique(data$subject))
  # Storage for the samples.
  # theta is the parameter values, mu is mean of normal distribution and
  # sigma2 is variance

  theta <- list(
    latent_mu = array(
      NA,
      dim = c(n_pars, n_subjects, 0),
      dimnames = list(parameters, NULL, NULL)
    ),
    param_mu = array(
      NA,
      dim = c(n_pars, 1, 0),
      dimnames = list(parameters, NULL, NULL)
    ),
    param_sigma2 <- array(
      NA,
      dim = c(n_pars, n_pars, 0),
      dimnames = list(parameters, parameters, NULL)
    )
  )
  # k_alpha from Algorithm 3, 2(b)
  k_half <- hyper$dof + n_pars - 1 + n_subjects
  # IG (inverse gamma) shape parameter, Algorithm 3, 2(c)
  v_shape <- (hyper$dof + n_pars) / 2
  # Sample the mixture variables' initial values.
  a_half <- 1 / stats::rgamma(n = n_pars, shape = 0.5, scale = 1)
  proposal <- list(
    means = array(dim = c(n_pars, n_subjects)),
    sigmas = array(dim = c(n_pars, n_pars, n_subjects))
  )
  sampler <- list(
    par_names = parameters,
    n_pars = n_pars,
    n_subjects = n_subjects,
    theta = theta,
    hyper = hyper,
    a_half = a_half,
    k_half = k_half,
    v_shape = v_shape
  )
  class(sampler) <- "pmwgs"
  sampler
}
