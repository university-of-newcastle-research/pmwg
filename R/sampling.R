#' Run a stage of the PMwG sampler
#'
#' Run one of burnin, adaptation of sampling phase from the PMwG
#' sampler. Each stage involves slightly different processes, so for the
#' full PMwG we need to run this three times.
#'
#' @param x A pmwgs object that has been initialised
#' @param stage The sampling stage to run. Must be one of 'burn', 'adapt' or
#'   'sample'. If not provided assumes that the stage should be 'burn'
#' @param iter The number of iterations to run for the sampler. For 'burn' and
#'   'sample' all iteration will run. However for 'adapt' if all subjects have
#'   enough unique samples to create the conditional distribution then it will
#'   finish early.
#' @param particles The default here is 1000 particles to be generated for each
#'   iteration, however during the sample phase this should be reduced.
#' @param display_progress Display a progress bar during sampling.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A pmwgs object with the newly generated samples in place.
#' @examples
#' # No example yet
#' @export
run_stage.pmwgs <- function(x, stage, iter = 1000, particles = 1000,  #nolint
                            display_progress = TRUE, ...) {
  # Test stage argument
  stage <- match.arg(stage, c("burn", "adapt", "sample"))
  if (stage == "sample") {
    eff <- try(create_efficient(x))
    if (class(eff) == "try-error") {
      cat("ERR01: An error was detected creating efficient dist\n")
      save.image("data/output/PMwG-error.RData")
      stop("Data saved in output directory under PMwG-error")
    }
    mix <- c(0.1, 0.2, 0.7)
    epsilon <- 1

  } else {
    mix <- c(0.5, 0.5, 0.0)
    epsilon <- 3
    eff <- list(prop_mean = NULL, prop_var = NULL)
  }
  # Test pmwgs object initialised
  try(if (is.null(x$init)) stop("pmwgs object has not been initialised"))

  # Display stage to screen
  msgs <- list(burn = "Phase 1: Burn in\n", adapt = "Phase 2: Adaptation\n",
               sample = "Phase 3: Sampling\n")
  cat(msgs[[stage]])

  # Build new sample storage
  stage_samples <- sample_store(x$par_names, x$n_subject, iters = iter)
  # create progress bar
  if (display_progress) {
    pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)
  }

  for (i in 1:iter) {

    if (display_progress) utils::setTxtProgressBar(pb, i)

    if (i == 1) store <- x$samples else store <- stage_samples
    pars <- new_group_pars(store, x)

    # Sample new particles for random effects.
    tmp <- lapply(
      X = 1:x$n_subjects,
      FUN = new_sample,
      data = x$data,
      num_particles = particles,
      parameters = pars,
      mix_ratio = mix,
      likelihood_func = x$ll_func,
      epsilon = epsilon,
      efficient_mu = eff$prop_mean,
      efficient_sig2 = eff$prop_var
    )
    sm <- array(unlist(tmp), dim = dim(pars$sm))

    # Store results.
    stage_samples$theta_mu[, i] <- pars$gm
    stage_samples$theta_sig[, , i] <- pars$gv
    stage_samples$last_theta_sig_inv <- pars$gvi
    stage_samples$alpha[, , i] <- sm
    stage_samples$idx <- i

    if (stage == "adapt") {
      if (check_adapted(stage_samples$alpha)) {
        print("Adapted - stopping early")
        break
      }
    }
  }
  if (display_progress) close(pb)
  update_sampler(x, stage_samples)
}


#' Generate a new particle
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
#' @param parameters A list containing:
#'          * the vector of means for the multivariate normal (gm),
#'          * A covariate matrix for the multivariate normal (gv)
#'          * An array of individual subject means (re proposals for latent
#'            variables) (sm)
#' @inheritParams numbers_from_ratio
#' @inheritParams check_efficient
#' @param likelihood_func A likelihood function for calculating log likelihood
#'   of samples
#' @param epsilon A scaling factor to reduce the variance on individual level
#'   parameter samples
#'
#' @return A single sample from the new proposals
#' @examples
#' # No example yet
#' @export
new_sample <- function(s, data, num_particles, parameters,
                       efficient_mu = NULL, efficient_sig2 = NULL,
                       mix_ratio = c(0.5, 0.5, 0.0),
                       likelihood_func = NULL,
                       epsilon = 1) {
  # Check for efficient proposal values if necessary
  check_efficient(mix_ratio, efficient_mu, efficient_sig2)
  e_mu <- efficient_mu[, s]
  e_sig2 <- efficient_sig2[, , s]
  mu <- parameters$gm
  sig2 <- parameters$gv
  subj_mu <- parameters$sm[, s]
  if (is.null(likelihood_func)) stop("likelihood_func is a required argument")

  # Create proposals for new particles
  proposals <- gen_particles(
    num_particles, mu, sig2, subj_mu,
    mix_ratio = mix_ratio,
    prop_mu = e_mu,
    prop_sig2 = e_sig2,
    epsilon = epsilon
  )
  # Put the current particle in slot 1.
  proposals[1, ] <- subj_mu

  # Density of data given random effects proposal.
  lw <- apply(
    proposals,
    1,
    likelihood_func,
    data = data[data$subject == s, ]
  )
  # Density of random effects proposal given population-level distribution.
  lp <- mvtnorm::dmvnorm(x = proposals, mean = mu, sigma = sig2, log = TRUE)
  # Density of proposals given proposal distribution.
  prop_density <- mvtnorm::dmvnorm(x = proposals, mean = subj_mu, sigma = sig2)
  # Density of efficient proposals
  if (mix_ratio[3] != 0) {
    eff_density <- mvtnorm::dmvnorm(
      x = proposals,
      mean = e_mu,
      sigma = e_sig2
    )
  } else {
    eff_density <- 0
  }

  lm <- log(mix_ratio[1] * exp(lp) +
    (mix_ratio[2] * prop_density) +
    (mix_ratio[3] * eff_density))
  # log of importance weights.
  l <- lw + lp - lm
  weights <- exp(l - max(l))
  proposals[sample(x = num_particles, size = 1, prob = weights), ]
}


#' Generate proposal particles.
#'
#' Generate particles for a particular subject from a mix of population level
#' (hierarchical) distribution, from the particles (containing individual level
#' distribution) and/or from the conditional on accepted individual level
#' particles, a more efficient prposal method.
#'
#' This function is used in burnin, adaptation and sampling using various
#' combinations of the arguments.
#'
#' @param mu A vector of means for the multivariate normal
#' @param sig2 A covariate matrix for the multivariate normal
#' @param particle A particle (re proposals for latent variables)
#' @inheritParams numbers_from_ratio
#' @param epsilon Reduce the variance for the individual level samples by this
#'   factor
#'
#' @return The new proposals
#' @examples
#' psamplers:::gen_particles(100, rep(0.2, 7), diag(rep(0.1, 7)), rep(0.3, 7))
#' @keywords internal
gen_particles <- function(num_particles,
                          mu,
                          sig2,
                          particle,
                          ...,
                          mix_ratio = c(0.5, 0.5, 0.0),
                          prop_mu = NULL,
                          prop_sig2 = NULL,
                          epsilon = 1) {
  particle_numbers <- numbers_from_ratio(mix_ratio, num_particles)
  # Generate proposal particles
  pop_particles <- particle_draws(particle_numbers[1], mu, sig2)
  ind_particles <- particle_draws(particle_numbers[2], particle, sig2 / epsilon)
  eff_particles <- particle_draws(particle_numbers[3], prop_mu, prop_sig2)
  particles <- rbind(pop_particles, ind_particles, eff_particles)
  colnames(particles) <- names(mu) # stripped otherwise.
  particles
}
