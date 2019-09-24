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
new_sample <- function(s, data, num_particles,
                       mu, sig2, particles, mix_ratio = 0.5) {
  # Create proposals for new particles
  proposals <- gen_particles(
    num_particles,
    mu,
    sig2,
    particles[, s],
    mix_ratio = mix_ratio
  )
  # Put the current particle in slot 1.
  proposals[1, ] <- particles[, s]
  # Density of data given random effects proposal.
  lw <- apply(
    proposals,
    1,
    lba_loglike,
    data = data[data$subject == s, ]
  )
  # Density of random effects proposal given population-level distribution.
  lp <- mvtnorm::dmvnorm(x = proposals, mean = mu, sigma = sig2, log = TRUE)
  # Density of proposals given proposal distribution.
  prop_density <- mvtnorm::dmvnorm(x = proposals,
                                   mean = particles[, s],
                                   sigma = sig2)
  lm <- log(mix_ratio * exp(lp) + (1 - mix_ratio) * prop_density)
  # log of importance weights.
  l <- lw + lp - lm
  weights <- exp(l - max(l))
  proposals[sample(x = num_particles, size = 1, prob = weights), ]
}


#' Generate proposal particles.
#'
#' Generate particles for a particular subject from a mix of population level (hierarchical)
#' distribution, from the particles (containing individual level distribution) and/or from
#' The conditional on accepted individual level particles, a more efficient prposal method.
#' This function is used in burnin, adaptation and sampling using various combinations of the arguments
#'
#' @param mu A vector of means for the multivariate normal
#' @param sig2 A covariate matrix for the multivariate normal
#' @param particle A particle (re proposals for latent variables)
#' @inheritParams numbers_from_ratio
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
                          proposal_means = NULL,
                          proposal_sigmas = NULL) {
  particle_numbers <- numbers_from_ratio(mix_ratio, num_particles)
  # Generate proposal particles
  population_particles <- particle_draws(particle_numbers[1], mu, sig2)
  randeffect_particles <- particle_draws(particle_numbers[2], particle, sig2)
  proposal_particles <- particle_draws(particle_numbers[3], proposal_means, proposal_sigmas)
  particles <- rbind(population_particles, randeffect_particles, proposal_particles)
  colnames(particles) <- names(mu) # stripped otherwise.
  particles
}
