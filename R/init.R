#' Initialise variables needed for individual loops within PMwG
#'
#' Takes the init list and further sets up individual loop variables
#'
#' @param init The list of parameter names to be used in the model
#' @param pmwg_args A list containing arguments to the PMwG sampling, including number of
#'   iterations for each stage.
#' @param prior A list containing priors for variance, means etc.
#' @param particles A multidimensional array containing accepted particlesfrom the sampling.
#'
#' @return A list of generated variables that can be modified after the fact
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
#' psamplers:::init_pmwg(c("b1", "b2", "b3", "A", "v1", "v2", "t0"), forstmann, args)
#' @export
gen_sample_pars <- function(init, pmwg_args, prior, particles) {
  # Sample population-level parameters.
  var_mu <- MASS::ginv(init$S * pts2_inv + prior$mu_sigma2_inv)
  mean_mu <- as.vector(
    var_mu %*% (pts2_inv %*% apply(particles, 1, sum))
  )
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  # New sample for mu.
  ptm <- mvtnorm::rmvnorm(
    1,
    mean_mu,
    chol_var_mu %*% t(chol_var_mu)
  )[1, ]
  names(ptm) <- init$par_names

  theta_temp <- particles - ptm
  cov_temp <- (theta_temp) %*% (t(theta_temp))
  B_half <- 2 * init$v_half * diag(1 / init$a_half) + cov_temp
  pts2 <- MCMCpack::riwish(init$k_half, B_half) # New sample for sigma.
  pts2_inv <- MASS::ginv(pts2)

  # Sample new mixing weights.
  init$a_half <- 1 / stats::rgamma(
    n = init$num_par,
    shape = init$v_shape,
    scale = 1 / (init$v_half + diag(pts2_inv) + init$A_half)
  )
  list(ptm = ptm, pts2 = pts2, pts2_inv = pts2_inv)
}
