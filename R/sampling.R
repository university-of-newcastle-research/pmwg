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
new_sample <- function(s, data, num_proposals,
                       mu, sig2, particles, mix_ratio = 0.5) {
  # Create proposals for new particles
  proposals <- psamplers::gen_proposals(
    num_proposals,
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
    psamplers::lba_loglike,
    data = data[data$subject == s, ]
  )
  # Density of random effects proposal given population-level distribution.
  lp <- mvtnorm::dmvnorm(x = proposals, mean = mu, sigma = sig2, log = TRUE)
  # Density of proposals given proposal distribution.
  prop_density <- mvtnorm::dmvnorm(x = proposals, mean = particles[, s], sigma = sig2)
  lm <- log(mix_ratio * exp(lp) + (1 - mix_ratio) * prop_density)
  # log of importance weights.
  l <- lw + lp - lm
  weights <- exp(l - max(l))
  proposals[sample(x = num_proposals, size = 1, prob = weights), ]
}

#' Generate proposal particles.
#'
#' Generate particles for a particular subject from a mix of either the
#' population level (hierarchical) distribution or from the particles
#' (containing individual level distribution).
#'
#' @param num_proposals A number representing the number of proposal particles to generate
#' @param mu A vector of means for the multivariate normal
#' @param sig2 A covariate matrix for the multivariate normal
#' @param particle A particle (re proposals for latent variables)
#' @param mix_ratio A float betwen 0 and 1 giving the ratio of particles to generate from the population level parameters vs the individual level parameters
#'
#' @return The new proposals
#' @examples
#' gen_proposals(100, rep(0.2, 7), diag(rep(0.1, 7)), rep(0.3, 7))
#' @keywords internal
gen_proposals <- function(num_proposals, mu, sig2, particle, mix_ratio = 0.5) {
  num_from_population <- rbinom(n = 1, size = num_proposals, prob = mix_ratio)
  if (num_from_population < 2) {
    num_from_population <- 2
  }
  if (num_from_population > (num_proposals - 2)) {
    num_from_population <- num_proposals - 2
  }
  num_from_subject <- num_proposals - num_from_population
  # Generate proposal particles
  population_proposals <- mvtnorm::rmvnorm(num_from_population, mu, sig2)
  subject_proposals <- mvtnorm::rmvnorm(num_from_subject, particle, sig2)
  proposals <- rbind(population_proposals, subject_proposals)
  colnames(proposals) <- names(mu) # stripped otherwise.
  proposals
}
