#' Run PMwG iterations
#'
#' Generate samples from either the initial proposal or from the last
#' iteration of the particles (according to a mixing value) and return
#' a weighted random sample from the new proposals
#' This uses the simplest, and slowest, proposals: a mixture of the
#' the population distribution and Gaussian around current random effect.
#'
#' @param s A number - the subject ID, also selects particles
#' @param data A data.frame containing variables for
#'        response time (\code{rt}), trial condition (\code{condition}),
#'        accuracy (\code{correct}) and subject (\code{subject}) which
#'        contains the data against which the particles are assessed
#' @param num_proposals A number representing the number of proposal particles to generate
#' @param mu A vector of means for the multivariate normal
#' @param sig2 A covariate matrix for the multivariate normal
#' @param particles An array of particles (re proposals for latent variables)
#' @param mix_ratio A float betwen 0 and 1 giving the ratio of particles to generate from the population level parameters vs the individual level parameters
#'
#' @return A single sample from the new proposals
#' @examples
#' # No example yet
#' @export
pmwg <- function(sample_data, pmwg_args, runner_args) {
  S <- sample_data$S
  data <- sample_data$data

  # Burnin
  # create progress bar
  pu <- runner_args$progress_update
  pb <- txtProgressBar(min = 0, max = pmwg_args$burnin_iter / pu, style = 3)

  for (i in 1:pmwg_args$sampling_iterations) {
    if (i %% pu == 0) {
      setTxtProgressBar(pb, i %/% pu)
    }
    # Sample population-level parameters.
    var_mu <- MASS::ginv(S * pts2_inv + prior_mu_sigma2_inv)
    mean_mu <- as.vector(var_mu %*% (pts2_inv %*% apply(particles, 1, sum)))
    chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
    ptm <- rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ] # New sample for mu.
    names(ptm) <- parameters

    theta_temp <- particles - ptm
    # cov.temp=array(0,dim=c(num_parameters,num_parameters))
    # for (j in 1:S) cov.temp=cov.temp+(theta.temp[,j])%*%t(theta.temp[,j])
    cov_temp <- (theta_temp) %*% (t(theta_temp))
    B_half <- 2 * v_half * diag(1 / a_half) + cov_temp
    pts2 <- riwish(k_half, B_half) # New sample for sigma.
    pts2_inv <- ginv(pts2)

    # Sample new mixing weights.
    a_half <- 1 / rgamma(n = num_parameters, shape = v_shape, scale = 1 / (v_half + diag(pts2_inv) + A_half))

    # Sample new particles for random effects.
    if (runner_args$cpus > 1) {
      tmp <- sfLapply(x = 1:S, fun = new_sample, data = data, num_proposals = num_particles, mu = ptm, sig2 = pts2, particles = particles)
    } else {
      tmp <- lapply(X = 1:S, FUN = new_sample, data = data, num_proposals = num_particles, mu = ptm, sig2 = pts2, particles = particles)
    }
    particles <- array(unlist(tmp), dim = dim(particles))

    # Store results.
    latent_theta_mu[, , i] <- particles
    param_theta_sigma2[, , i] <- pts2
    param_theta_mu[, i] <- ptm
  }
  close(pb)
  if (runner_args$cpus > 1) sfStop()
}
