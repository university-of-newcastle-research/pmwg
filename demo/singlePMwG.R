library(mvtnorm) ## for the multivariate normal.
library(MASS) ## for matrix inverse.
library(MCMCpack) ## for the inverse wishart random numbers.
library(psamplers)
rm(list = ls())

# Vars used for controlling the run
sampler <- pmwgs(forstmann, c("b1", "b2", "b3", "A", "v1", "v2", "t0"))


start_points_mu <- c(.2, .2, .2, .4, .3, 1.3, -2)
start_points_sig2 <- diag(rep(.01, init$num_par))

# Storing the proposal values for the conditional Monte Carlo

# Make single-iteration-sized versions, for easier reading of code below.
ptm <- init$param_theta_mu[, 1]
pts2 <- init$param_theta_sigma2[, , 1] # nolint

# Start points for the population-level parameters only. Hard coded here, just
# for convenience.
# Weird subscripts maintains the naming.
# GC: Check the start points are correct?
ptm[1:init$num_par] <- start_points_mu
pts2 <- start_points_sig2 # Who knows??
# Because this is calculated near the end of the main loop, needs initialising for iter=1.
pts2_inv <- ginv(pts2)

# Priors.
prior <- list(
  mu_mean = rep(0, init$num_par),
  mu_sigma2 = diag(rep(1, init$num_par))
)
# Things I save rather than re-compute inside the loops.
prior$mu_sigma2_inv <- ginv(prior$mu_sigma2)

# Sample the initial values for the random effects. Algorithm is same
# as for the main resampling down below.
particles <- array(dim = c(length(ptm), init$S))
num_particles <- pmwg_args$burn_particles
cat("Sampling Initial values for random effects\n")
pb <- txtProgressBar(min = 0, max = init$S, style = 3)
for (s in 1:init$S) {
  setTxtProgressBar(pb, s)
  proposals <- rmvnorm(num_particles, ptm, pts2)
  lw <- apply(
    proposals,
    1,
    pmwg_args$likelihood_func,
    data = data[data$subject == s, ]
  )
  weight <- exp(lw - max(lw))
  particles[, s] <- proposals[
    sample(x = num_particles, size = 1, prob = weight),
  ]
}
close(pb)

cat("Phase 1: Burn in\n")
# create progress bar
pb <- txtProgressBar(
  min = 0,
  max = pmwg_args$burn_iter / progress_update,
  style = 3
)

for (i in 1:pmwg_args$burn_iter) {
  if (i %% progress_update == 0) {
    setTxtProgressBar(pb, i %/% progress_update)
  }

  single_iter <- gen_sample_pars(init, pmwg_args, prior, particles)
  ptm <- single_iter$ptm
  pts2 <- single_iter$pts2

  # Sample new particles for random effects.
  tmp <- lapply(
    X = 1:init$S,
    FUN = new_sample,
    data = data,
    num_particles = num_particles,
    mu = ptm,
    sig2 = pts2,
    particles = particles,
    mix_ratio = c(0.5, 0.5, 0.0)
  )
  particles <- array(unlist(tmp), dim = dim(particles))

  # Store results.
  init$latent_theta_mu[, , i] <- particles # nolint
  init$param_theta_sigma2[, , i] <- pts2 # nolint
  init$param_theta_mu[, i] <- ptm
}
close(pb)

# Adaptation Phase
num_particles <- pmwg_args$adapt_particles
cat("Phase 2: Adaptation\n")
# create progress bar
pb <- txtProgressBar(
  min = 0,
  max = pmwg_args$adapt_maxiter / progress_update,
  style = 3
)
pmwg_args$adapted <- FALSE
pmwg_args$adapt_iter <- pmwg_args$burn_iter + pmwg_args$adapt_maxiter

for (i in 1:pmwg_args$adapt_maxiter) {
  if (i %% progress_update == 0) {
    setTxtProgressBar(pb, i %/% progress_update)
  }

  single_iter <- gen_sample_pars(init, pmwg_args, prior, particles)
  ptm <- single_iter$ptm
  pts2 <- single_iter$pts2

  # Sample new particles for random effects.
  tmp <- lapply(
    X = 1:init$S,
    FUN = new_sample,
    data = data,
    num_particles = num_particles,
    mu = ptm,
    sig2 = pts2,
    particles = particles,
    mix_ratio = c(0.5, 0.5, 0.0)
  )
  particles <- array(unlist(tmp), dim = dim(particles))

  # Store results.
  init$latent_theta_mu[, , pmwg_args$burn_iter + i] <- particles # nolint
  init$param_theta_sigma2[, , pmwg_args$burn_iter + i] <- pts2 # nolint
  init$param_theta_mu[, pmwg_args$burn_iter + i] <- ptm

  # Check adaptive phase ended
  # Check the number of unique 'A' values since the end of the burnin
  # Get the length of each vector of unique values
  # Test the if every length is greater than 20
  pmwg_adapted <- all(
    lapply(
      apply(
        init$latent_theta_mu["A", , pmwg_args$burn_iter:(pmwg_args$burn_iter + i)],
        1,
        unique
      ),
      length
    ) > 20
  )
  if (pmwg_adapted) {
    pmwg_args$adapted <- TRUE
    pmwg_args$adapt_iter <- pmwg_args$burn_iter + i
    break
  }
}
close(pb)

save.image("data/output/PMwG.RData")
#Create conditional means/variances
for (s in 1:init$S) {
  cparms <- conditional_parms(init,
                              ptm,
                              pts2,
                              s,
                              pmwg_args$burn_iter + 1,
                              pmwg_args$burn_iter + pmwg_args$adapt_iter)
  proposal_means[, s] <- cparms$cmeans
  proposal_sigmas[, , s] <- cparms$cvars
}

# Sampling Phase
num_particles <- pmwg_args$sample_particles
cat("Phase 3: Sampling\n")
# create progress bar
pb <- txtProgressBar(
  min = 0,
  max = pmwg_args$sample_iter / progress_update,
  style = 3
)

for (i in 1:pmwg_args$sample_iter) {
  if (i %% progress_update == 0) {
    setTxtProgressBar(pb, i %/% progress_update)
  }

  single_iter <- gen_sample_pars(init, pmwg_args, prior, particles)
  ptm <- single_iter$ptm
  pts2 <- single_iter$pts2

  # Sample new particles for random effects.
  tmp <- lapply(
    X = 1:init$S,
    FUN = new_sample,
    data = data,
    num_particles = num_particles,
    mu = ptm,
    sig2 = pts2,
    particles = particles,
    mix_ratio = c(0.1, 0.2, 0.7)
  )
  particles <- array(unlist(tmp), dim = dim(particles))

  # Store results.
  init$latent_theta_mu[, , pmwg_args$adapt_iter + i] <- particles # nolint
  init$param_theta_sigma2[, , pmwg_args$adapt_iter + i] <- pts2 # nolint
  init$param_theta_mu[, pmwg_args$adapt_iter + i] <- ptm
}
close(pb)

save.image("data/output/PMwG.RData")
