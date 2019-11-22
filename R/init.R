#' Initialise variables needed for PMwG with default values
#'
#' Takes the initial parameter set and creates default values for a
#' range of PMwG necessary variables.
#'
#' @param parameters The list of parameter names to be used in the model
#' @param data The data.frame containing empirical data to be modelled. Assumed to contain
#'   at least one column called subject whose elements are unique identifiers for each subject.
#' @param pmwg_args A list containing arguments to the PMwG sampling, including number of
#'   iterations for each stage.
#'
#' @return A list of generated variables that can be modified after the fact
#' @examples
#' #Load Forstmann et al.'s data.
#' args <- list(
#'   "burn_particles" = 1000,
#'   "burn_iter" = 500,
#'   "adapt_particles" = 100,
#'   "adapt_maxiter" = 5000,
#'   "sample_particles" = 100,
#'   "sample_iter" = 1000,
#'   "likelihood_func" = lba_loglike
#' )
#' psamplers:::init_pmwg(c("b1", "b2", "b3", "A", "v1", "v2", "t0"), forstmann, args)
#' @export
init_pmwg <- function(parameters, data, pmwg_args) {
  init <- list()
  init$par_names <- parameters
  init$num_par <- length(parameters)
  # Tuning settings for the Gibbs steps
  init$v_half <- 2
  init$A_half <- 1
  init$S <- length(unique(data$subject))
  # Storage for the samples.
  # theta is the parameter values, mu is mean of normal distribution and
  # sigma2 is variance
  max_iter <- sum(pmwg_args$burn_iter,
                  pmwg_args$adapt_maxiter,
                  pmwg_args$sample_iter)
  init$latent_theta_mu <- array(
    NA,
    dim = c(init$num_par, init$S, max_iter),
    dimnames = list(parameters, NULL, NULL)
  )
  init$param_theta_mu <- init$latent_theta_mu[, 1, ]
  init$param_theta_sigma2 <- array(
    NA,
    dim = c(init$num_par, init$num_par, pmwg_args$sample_iter),
    dimnames = list(parameters, parameters, NULL)
  )
  init$k_half <- init$v_half + init$num_par - 1 + init$S
  init$v_shape <- (init$v_half + init$num_par) / 2
  # Sample the mixture variables' initial values.
  init$a_half <- 1 / rgamma(n = init$num_par, shape = 0.5, scale = 1)
	init
}
