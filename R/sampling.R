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
#' @param mu A vector of means for the multivariate normal
#' @param sig2 A covariate matrix for the multivariate normal
#' @param particles An array of particles (re proposals for latent variables)
#' @inheritParams numbers_from_ratio
#'
#' @return A single sample from the new proposals
#' @examples
#' # No example yet
#' @export
new_sample <- function(s, data, num_particles,
                       mu, sig2, particles,
                       efficient_mu = NULL, efficient_sig2 = NULL,
                       mix_ratio = c(0.5, 0.5, 0.0)) {
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
  #      %computing the log density of the LBA given the particles of random
  #    %effects
  #
  #    logw_first=sum(lw_reshape);
  #
  #    %computing the log of p(\alpha|\theta) and density of the proposal for
  #    %burn in and initial sampling stage (count<=switch_num) and sampling
  #    %stage (count>switch_num)
  #
  #    logw_second=(logmvnpdf(rnorm_theta,param.theta_mu,chol_covmat*chol_covmat'));
  #    if  sum(count>switch_num)==num_subjects
  #        logw_third=log(w_mix.*mvnpdf(rnorm_theta,cond_mean',chol_cond_var*chol_cond_var')+...
  #            (1-w_mix).*mvnpdf(rnorm_theta,param.theta_mu,chol_covmat*chol_covmat'));
  #    else
  #        logw_third=log(w_mix.*mvnpdf(rnorm_theta,reference_par,(epsilon^2).*(chol_covmat*chol_covmat'))+...
  #            (1-w_mix).*mvnpdf(rnorm_theta,param.theta_mu,chol_covmat*chol_covmat'));
  #    end
  #    logw=logw_first'+logw_second'-logw_third;

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
  prop_density <- mvtnorm::dmvnorm(
    x = proposals,
    mean = particles[, s],
    sigma = sig2
  )
  # Density of efficient proposals
  eff_density <- mvtnorm::dmvnorm(x = proposals,
                                  mean = efficient_mu,
                                  sigma = efficient_sig2
                                 )

  lm <- log(mix_ratio[1] * exp(lp) +
        (mix_ratio[2] * prop_density) +
        (mix_ratio[3] * le))
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
