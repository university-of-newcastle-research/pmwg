#' Initialise a PMwG sampler and return a pmwgs object
#'
#' Takes the initial parameter set and creates default values for a
#' range of PMwG necessary variables.
#'
#' @param data The data.frame containing empirical data to be modelled. Assumed
#'   to contain at least one column called subject whose elements are unique
#'   identifiers for each subject.
#' @param parameters The list of parameter names to be used in the model
#' @param pmwg_args A list containing arguments to the PMwG sampling, including
#'   number of iterations for each stage.
#'
#' @return A pmwgs object that can be run, plotted and more
#' @examples
#' # Load Forstmann et al.'s data.
#' args <- list(
#'   "burn_particles" = 1000,
#'   "burn_iter" = 500,
#'   "adapt_particles" = 100,
#'   "adapt_maxiter" = 5000,
#'   "sample_particles" = 100,
#'   "sample_iter" = 1000,
#'   "likelihood_func" = lba_loglike
#' )
#' pmwgs(c("b1", "b2", "b3", "A", "v1", "v2", "t0"), forstmann)
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
  max_iter <- burn[2] + adapt[2] + sample[2]

  latent_theta_mu <- array(
    NA,
    dim = c(n_pars, n_subjects, max_iter),
    dimnames = list(parameters, NULL, NULL)
  )
  param_theta_mu <- latent_theta_mu[, 1, ]
  param_theta_sigma2 <- array(
    NA,
    dim = c(n_pars, n_pars, max_iter),
    dimnames = list(parameters, parameters, NULL)
  )
  k_half <- hyper$dof + n_pars - 1 + n_subjects
  v_shape <- (hyper$dof + n_pars) / 2
  # Sample the mixture variables' initial values.
  a_half <- 1 / stats::rgamma(n = n_pars, shape = 0.5, scale = 1)
  proposal <- list(
    means = array(dim = c(n_pars, n_subjects)),
    sigmas = array(dim = c(n_pars, n_pars, n_subjects))
  )
  sampler <- list(
    par_names = parameters,
    num_par = num_par,
    n_subjects = n_subjects,
    param_theta_mu = param_theta_mu,
    param_theta_sigma2 = param_theta_sigma2,
    latent_theta_mu = latent_theta_mu,
    hyper = hyper,
    a_half = a_half,
    k_half = k_half,
    v_shape = v_shape
  )
  class(sampler) <- "pmwgs"
  sampler
}
