#' Paralellisable particle sampling function (cond/orig/real/log).
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
#'
#' @return A single sample from the new proposals
#' @examples
#' # No example yet
#' @export
new_sample <- function(s, data, num_proposals,
                       mu, sig2, particles, mix_ratio = 0.5) {
  # Choose the number of particles to generate from one of two methods
  # Restrict the number of each type of proposal to be at least 2 to avoid
  # degenerate arrays
  num_from_population <- rbinom(n = 1, size = num_proposals, prob = mix_ratio)
  if (num_from_population < 2) {
    num_from_population <- 2
  }
  if (num_from_population > (num_proposals - 2)) {
    num_from_population <- num_proposals - 2
  }
  num_from_individual <- num_proposals - num_from_population
  # Generate proposal particles
  population_proposals <- mvtnorm::rmvnorm(num_from_population, mu, sig2)
  # Generate other proposals using the individual particles for this subject
  individual_proposals <- mvtnorm::rmvnorm(
    num_from_individual,
    particles[, s],
    sig2
  )
  proposals <- rbind(population_proposals, individual_proposals)
  colnames(proposals) <- names(mu) # stripped otherwise.
  proposals[1, ] <- particles[, s] # Put the current particle in slot 1.
  # Density of data given random effects proposal.
  lw <- apply(
    proposals,
    1,
    psamplers::lba_loglike,
    data = data[data$subject == s, ]
  )
  # Density of random effects proposal given population-level distribution.
  lp <- dmvnorm(x = proposals, mean = mu, sigma = sig2, log = TRUE)
  # Density of proposals given proposal distribution.
  lm <- log(mix_ratio * exp(lp) + (1 - mix_ratio) * dmvnorm(x = proposals, mean = particles[, s], sigma = sig2))
  l <- lw + lp - lm # log of importance weights.
  weight <- exp(l - max(l))
  proposals[sample(x = num_proposals, size = 1, prob = weight), ]
}
