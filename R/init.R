#' Initialise values for the random effects
#'
#' Initialise the random effects for each subject using MCMC.
#'
#' Before sampling can start the Particle Metropolis within Gibbs sampler needs
#' initial values for the random effects. The \code{init} function generates
#' these values using a Monte Carlo algorithm. One alternative methods would be
#' setting the initial values randomly.
#'
#' Optionally takes starting values for the model parameters and the variance /
#' covariance matrix. All arrays must match the appropriate shape.
#'
#' For example, with 5 parameters and 10 subjects, the model parameter start
#' means must be a vector of length 5 and the covariance matrix must be an array
#' of 5 x 5.
#'
#' If the start_mu and start_sig arguments are left at the default (NULL) then
#' start_mu will be sampled from a normal distribution with mean as the prior
#' mean for eac variable and sd as the square of the variance from the prior
#' covariance matrix. start_sig by default is sampled from an inverse wishart
#' (IW) distribution. For a model with the number of parameters N the degrees of
#' freedom of the IW distribution is set to N*3 and the scale matrix is the
#' identity matrix of size NxN.
#'
#' @param pmwgs The sampler object that provides the parameters.
#' @param start_mu An array of starting values for the group means
#' @param start_sig An array of starting values for the group covariance matrix
#' @param display_progress Display a progress bar during sampling
#' @param particles The number of particles to generate in initialisation
#'
#' @return The sampler object but with initial values set for \code{theta_mu},
#'   \code{theta_sig}, \code{alpha} and other values for the first sample.
#' @examples
#' lba_ll <- function(x, data) {
#'   x <- exp(x)
#'   if (any(data$rt < x["t0"])) {
#'     return(-1e10)
#'   }
#'   sum(
#'     log(
#'       rtdists::dLBA(
#'         rt = data$rt,
#'         response = data$correct,
#'         A = x["A"],
#'         b = x["A"] + x[c("b1", "b2", "b3")][data$condition],
#'         t0 = x["t0"],
#'         mean_v = x[c("v1", "v2")],
#'         sd_v = c(1, 1),
#'         silent = TRUE
#'       )
#'     )
#'   )
#' }
#' sampler <- pmwgs(
#'   forstmann,
#'   c("b1", "b2", "b3", "A", "v1", "v2", "t0"),
#'   lba_ll
#' )
#' sampler <- init(sampler, particles=10)
#' @export
init <- function(pmwgs, start_mu = NULL, start_sig = NULL,
                 display_progress = TRUE, particles = 100) {
  if (is.null(attr(pmwgs, "class"))) {
    message("ERROR: No object to add start points to")
  }
  # If no starting point for group mean just sample from prior
  if (is.null(start_mu)) {
    start_mu <- stats::rnorm(
      pmwgs$n_pars,
      mean = pmwgs$prior$theta_mu_mean,
      sd = sqrt(diag(pmwgs$prior$theta_mu_var)))
  }
  # If no starting point for group var just sample from inverse wishart
  if (is.null(start_sig)) {
    start_sig <- riwish(pmwgs$n_pars * 3, diag(pmwgs$n_pars))
  }
  n_particles <- particles
  # Sample the mixture variables' initial values.
  a_half <- 1 / stats::rgamma(n = pmwgs$n_pars, shape = 0.5, scale = 1)
  # Create and fill initial random effects for each subject
  alpha <- array(NA, dim = c(pmwgs$n_pars, pmwgs$n_subjects))
  if (display_progress) {
    message("MESSAGE: Sampling Initial values for random effects")
    pb <- utils::txtProgressBar(min = 0, max = pmwgs$n_subjects, style = 3)
  }
  likelihoods <- array(NA_real_, dim = c(pmwgs$n_subjects))
  for (s in 1:pmwgs$n_subjects) {
    if (display_progress) utils::setTxtProgressBar(pb, s)
    particles <- mvtnorm::rmvnorm(n_particles, start_mu, start_sig)
    colnames(particles) <- rownames(pmwgs$samples$theta_mu) # preserve par names
    lw <- apply(
      particles,
      1,
      pmwgs$ll_func,
      data = pmwgs$data[pmwgs$data$subject == pmwgs$subjects[s], ]
    )
    weight <- exp(lw - max(lw))
    idx <- sample(x = n_particles, size = 1, prob = weight)
    alpha[, s] <- particles[idx, ]
    likelihoods[s] <- lw[idx]
  }
  if (display_progress) close(pb)
  pmwgs$init <- TRUE
  pmwgs$samples$theta_mu[, 1] <- start_mu
  pmwgs$samples$theta_sig[, , 1] <- start_sig
  pmwgs$samples$alpha[, , 1] <- alpha
  pmwgs$samples$last_theta_sig_inv <- MASS::ginv(start_sig)
  pmwgs$samples$subj_ll[, 1] <- likelihoods
  pmwgs$samples$a_half[, 1] <- a_half
  pmwgs$samples$idx <- 1
  pmwgs
}


#' Gibbs step of the Particle Metropolis within Gibbs sampler
#'
#' Samples new \code{theta_mu} and \code{theta_sig} using Gibbs sampling
#'
#' @param sampler The pmwgs object from which to generate the new group
#'   parameters.
#'
#' @return A new sample for \code{theta_mu}, \code{theta_sig} and some new
#'   mixing weights in a list for use in the Particle Metropolis step.
#' @keywords internal
gibbs_step <- function(sampler) {
  # Get single iter versions, tmu = theta_mu, tsig = theta_sig
  last <- last_sample(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior

  # Here mu is group mean, so we are getting mean and variance
  var_mu <- MASS::ginv(
    sampler$n_subjects * last$tsinv + prior$theta_mu_invar
  )
  mean_mu <- as.vector(
    var_mu %*% (last$tsinv %*% apply(last$alpha, 1, sum) +
                prior$theta_mu_invar %*% prior$theta_mu_mean)
  )
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  # New sample for mu.
  tmu <- mvtnorm::rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
  names(tmu) <- sampler$par_names

  # New values for group var
  theta_temp <- last$alpha - tmu
  cov_temp <- (theta_temp) %*% (t(theta_temp))
  B_half <- 2 * hyper$v_half * diag(1 / last$a_half) + cov_temp # nolint
  tsig <- riwish(hyper$k_half, B_half) # New sample for group variance
  tsinv <- MASS::ginv(tsig)

  # Sample new mixing weights.
  a_half <- 1 / stats::rgamma(
    n = sampler$n_pars,
    shape = hyper$v_shape,
    scale = 1 / (hyper$v_half * diag(tsinv) + hyper$A_half)
  )
  list(
    tmu = tmu,
    tsig = tsig,
    tsinv = tsinv,
    a_half = a_half,
    alpha = last$alpha
  )
}
