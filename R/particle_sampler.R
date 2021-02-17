#' Create a PMwG sampler and return the created object
#'
#' This function takes a few necessary elements for creating a PMwG sampler.
#' Each pmwgs object is required to have a dataset, a list of parameter names,
#' a log likelihood function and optionally a prior for the model parameters.
#'
#' @param data A data frame containing empirical data to be modelled. Assumed
#'   to contain at least one column called subject whose elements are unique
#'   identifiers for each subject. Can be any of \code{data.frame},
#'   \code{data.table} or \code{tibble}, or any other data frame like object
#'   that can have subsets created in an identical way.
#' @param pars The list of parameter names to be used in the model
#' @param ll_func A log likelihood function that given a list of parameter
#'   values and a data frame (or other data store) containing subject data will
#'   return the log likelihood of \code{x} given \code{data}.
#' @param prior Specification of the prior distribution for the model
#'   parameters. It should be a list with two elements, \code{theta_mu_mean} and
#'   \code{theta_mu_var} which fully specify the prior distribution. If left as
#'   the default (NULL) then the \code{theta_mu_mean} will be all zeroes and
#'   \code{theta_mu_var} will be 1 on the diagonal and zero elsewhere.
#'
#' @return A pmwgs object that is ready to be initialised and sampled.
#' @example examples/pmwgs.R
#' @export
pmwgs <- function(data, pars, ll_func, prior = NULL) {
  # Descriptives
  n_pars <- length(pars)
  if (!"subject" %in% colnames(data)) {
    stop("Data must have a column named 'subject'")
  }
  subjects <- unique(data$subject)
  n_subjects <- length(subjects)
  # Tuning settings for the Gibbs steps
  # Hyperparameters
  v_half <- 2 # hyperparameter on Σ prior (Half-t degrees of freedom)
  A_half <- 1 # hyperparameter on Σ prior (Half-t scale) #nolint
  # k_alpha from Algorithm 3, 2(b)
  k_half <- v_half + n_pars - 1 + n_subjects
  # Inverse Gamma shape parameter, Algorithm 3, 2(c)
  v_shape <- (v_half + n_pars) / 2
  # Storage for the samples.
  samples <- sample_store(pars, subjects)
  # Checking and default priors
  prior_default <- list(
    theta_mu_mean = rep(0, n_pars),
    theta_mu_var = diag(rep(1, n_pars))
  )
  if (is.null(prior)) {
    prior <- prior_default
  }
  else {
    if (!identical(names(prior), names(prior_default))) {
      stop(paste(
        "pmwgs prior should be a list with two elements,",
        "`theta_mu_mean`, a vector that is the prior for the mean of the",
        "group-level mean parameters and `theta_mu_var`, a covariance matrix",
        "that is the prior for the variance of the group-level mean",
        "parameters"))
    }
    if (!identical(lapply(prior, dim), lapply(prior_default, dim)) |
        !identical(lapply(prior, length), lapply(prior_default, length))) {
      stop(paste(
        "pmwgs prior list elements specified with incorrect shape.",
        "`theta_mu_mean` should have length equal to the pars argument.",
        "`theta_mu_var` should have dim equal to N x N where N is length of",
        "pars"))
    }
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- MASS::ginv(prior$theta_mu_var)

  sampler <- list(
    data = data,
    par_names = pars,
    n_pars = n_pars,
    n_subjects = n_subjects,
    subjects = subjects,
    prior = prior,
    ll_func = ll_func,
    samples = samples,
    init = FALSE
  )
  attr(sampler, "v_half") <- v_half
  attr(sampler, "A_half") <- A_half
  attr(sampler, "k_half") <- k_half
  attr(sampler, "v_shape") <- v_shape
  class(sampler) <- "pmwgs"
  sampler
}

#' Test whether object is a pmwgs
#'
#' Checks whether object is a Particle Metropolis with Gibbs sampler
#'
#' @param x An object to test
#'
#' @return logical, whether object inherits from pmwgs
#' @examples
#' if (is.pmwgs(sampled_forstmann)) {
#'   print("sampled_forstmann object is a pmwgs")
#' }
#' @export
is.pmwgs <- function(x) inherits(x, "pmwgs")  # nolint

stages <- c("init", "burn", "adapt", "sample")
