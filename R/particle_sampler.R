#' Initialise a PMwG sampler and return a pmwgs object
#'
#' Takes the initial parameter set and creates default values for a
#' range of PMwG necessary variables.
#'
#' @param data The data.frame containing empirical data to be modelled. Assumed
#'   to contain at least one column called subject whose elements are unique
#'   identifiers for each subject.
#' @param parameters The list of parameter names to be used in the model
#' @param llfunc A log likelihood function that given a list of parameter values
#'   and a data.frame (or other data store) containing subject data will return
#'   the log likelihood of \code{x} given \code{data}.
#' @param prior Specification of the prior distribution for mu and sigma2
#'
#' @return A pmwgs object that can be run, plotted and more
#' @example examples/pmwgs.R
#' @export
pmwgs <- function(data, parameters, llfunc, prior = NULL) {
  n_pars <- length(parameters)
  # Tuning settings for the Gibbs steps
  hyper <- list(
    dof = 2,  # hyperparameter on prior (Half-t degrees of freedom)
    scale = 1 # hyperparameter on prior (Half-t scale)
  )
  n_subjects <- length(unique(data$subject))
  # Storage for the samples.
  # theta is the parameter values, mu is mean of normal distribution and
  # sigma2 is variance

  # Generate a list of sample storage arrays of size 1 (for start points)
  samples <- sample_store(parameters, n_subjects, 1)
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
  if (is.null(prior)) {
    prior <- list(group_mean = rep(0, n_pars), group_var = diag(rep(1, n_pars)))
  }
  # Things I save rather than re-compute inside the loops.
  prior$group_var_inv <- MASS::ginv(prior$group_var)

  sampler <- list(
    data = data,
    par_names = parameters,
    n_pars = n_pars,
    n_subjects = n_subjects,
    prior = prior,
    llfunc = llfunc,
    samples = samples,
    hyper = hyper,
    a_half = a_half,
    k_half = k_half,
    v_shape = v_shape
  )
  class(sampler) <- "pmwgs"
  sampler
}
