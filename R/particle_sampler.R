#' Initialise a PMwG sampler and return a pmwgs object
#'
#' Takes the initial parameter set and creates default values for a
#' range of PMwG necessary variables.
#'
#' @param data The data.frame containing empirical data to be modelled. Assumed
#'   to contain at least one column called subject whose elements are unique
#'   identifiers for each subject.
#' @param pars The list of parameter names to be used in the model
#' @param ll_func A log likelihood function that given a list of parameter
#'  values and a data.frame (or other data store) containing subject data will
#'  return the log likelihood of \code{x} given \code{data}.
#' @param prior Specification of the prior distribution for mu and sigma2
#' @param ... Other tuning parameters that you want to specify for the sampler.
#'
#' @return A pmwgs object that can be run, plotted and more
#' @example examples/pmwgs.R
#' @export
pmwgs <- function(data, pars, ll_func, prior = NULL, ...) {
  # Descriptives
  n_pars <- length(pars)
  n_subjects <- length(unique(data$subject))
  # Tuning settings for the Gibbs steps
  # Hyperparameters
  v_half <- 2 # hyperparameter on Σ prior (Half-t degrees of freedom)
  A_half <- 1 # hyperparameter on Σ prior (Half-t scale) #nolint
  # k_alpha from Algorithm 3, 2(b)
  k_half <- v_half + n_pars - 1 + n_subjects
  # IG (inverse gamma) shape parameter, Algorithm 3, 2(c)
  v_shape <- (v_half + n_pars) / 2
  # Sample the mixture variables' initial values.
  a_half <- 1 / stats::rgamma(n = n_pars, shape = 0.5, scale = 1)
  # Storage for the samples.
  # theta is the parameter values, mu is mean of normal distribution and
  # sigma2 is variance
  # Generate a list of sample storage arrays of size 1 (for start points)
  samples <- sample_store(pars, n_subjects, 1)
  proposal <- list(
    means = array(dim = c(n_pars, n_subjects)),
    sigmas = array(dim = c(n_pars, n_pars, n_subjects))
  )
  if (is.null(prior)) {
    prior <- list(theta_mu = rep(0, n_pars), theta_sig = diag(rep(1, n_pars)))
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_sig_inv <- MASS::ginv(prior$theta_sig)

  sampler <- list(
    data = data,
    par_names = pars,
    n_pars = n_pars,
    n_subjects = n_subjects,
    prior = prior,
    ll_func = ll_func,
    samples = samples
  )
  attr(sampler, "v_half") <- v_half
  attr(sampler, "A_half") <- A_half
  attr(sampler, "a_half") <- a_half
  attr(sampler, "k_half") <- k_half
  attr(sampler, "v_shape") <- v_shape
  class(sampler) <- "pmwgs"
  sampler
}
