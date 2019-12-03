library(mvtnorm) ## for the multivariate normal.
library(MASS) ## for matrix inverse.
library(MCMCpack) ## for the inverse wishart random numbers.
library(psamplers)
rm(list = ls())

# Vars used for controlling the run
restart <- FALSE
cpus <- 4
restart_file <- "data/output/restart.RData"
progress_update <- 10

pmwg_args <- list(
  "burn_particles" = 1000,
  "burn_iter" = 500,
  "adapt_particles" = 100,
  "adapt_maxiter" = 5000,
  "sample_particles" = 100,
  "sample_iter" = 1000,
  "likelihood_func" = lba_loglike
)

# Load Forstmann et al.'s data.
data <- forstmann

init <- init_pmwg(c("b1", "b2", "b3", "A", "v1", "v2", "t0"), data, pmwg_args)
start_points_mu <- c(.2, .2, .2, .4, .3, 1.3, -2)
start_points_sig2 <- diag(rep(.01, init$num_par))

# Make single-iteration-sized versions, for easier reading of code below.
ptm <- init$param_theta_mu[, 1]
pts2 <- init$param_theta_sigma2[, , 1] # nolint

# Start points for the population-level parameters only. Hard coded here, just
# for convenience.
# Weird subscripts maintains the naming.
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

if (restart) {
  cat("\nRestarting from saved run.\n")
  load(restart_file)
  # Check a couple of things.
  if ((dim(particles)[1] != init$num_par) || (dim(particles)[2] != init$S)) { # nolint
    stop("Restart does not match size (subjects).")
  }
}

if (cpus > 1) {
  # nolint start
  # You need the suggested package for this function
  library(snowfall)

  sfInit(parallel = TRUE, cpus = cpus)
  sfClusterSetupRNG()
  sfLibrary(rtdists)
  sfLibrary(mvtnorm)
  sfLibrary(MASS) # For matrix inverse.
  sfLibrary(MCMCpack) # For the inverse Wishart random numbers.

  sfExportAll(except = c(
    "init"
  ))
  # nolint end
}

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

  # Sample new particles for random effects.
  if (cpus > 1) {
    tmp <- sfLapply( # nolint
      x = 1:init$S,
      fun = new_sample,
      data = data,
      num_particles = num_particles,
      mu = single_iter$ptm,
      sig2 = single_iter$pts2,
      particles = particles,
      mix_ratio = c(0.5, 0.5, 0.0)
    )
  } else {
    tmp <- lapply(
      X = 1:init$S,
      FUN = new_sample,
      data = data,
      num_particles = num_particles,
      mu = single_iter$ptm,
      sig2 = single_iter$pts2,
      particles = particles,
      mix_ratio = c(0.5, 0.5, 0.0)
    )
  }
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
  max = pmwg_args$burn_iter / progress_update,
  style = 3
)

for (i in 1:pmwg_args$adapt_maxiter) {
  if (i %% progress_update == 0) {
    setTxtProgressBar(pb, i %/% progress_update)
  }

  single_iter <- gen_sample_pars(init, pmwg_args, prior, particles)

  # Sample new particles for random effects.
  if (cpus > 1) {
    tmp <- sfLapply( # nolint
      x = 1:init$S,
      fun = new_sample,
      data = data,
      num_particles = num_particles,
      mu = single_iter$ptm,
      sig2 = single_iter$pts2,
      particles = particles,
      mix_ratio = c(0.5, 0.5, 0.0)
    )
  } else {
    tmp <- lapply(
      X = 1:init$S,
      FUN = new_sample,
      data = data,
      num_particles = num_particles,
      mu = single_iter$ptm,
      sig2 = single_iter$pts2,
      particles = particles,
      mix_ratio = c(0.5, 0.5, 0.0)
    )
  }
  particles <- array(unlist(tmp), dim = dim(particles))

  # Store results.
  init$latent_theta_mu[, , pmwg_args$burn_iter + i] <- particles # nolint
  init$param_theta_sigma2[, , pmwg_args$burn_iter + i] <- pts2 # nolint
  init$param_theta_mu[, pmwg_args$burn_iter + i] <- ptm

  # Check adaptive phase ended
  # Check the number of unique 'A' values since the end of the burnin (apply call)
  # Get the length of each vector of uniqie values
  # Test the if every length is greater than 20
  pmwg_adapted <- all(
    lapply(
      apply(
        init$latent_theta_mu["A", , pmwg_args$burn_iter:pmwg_args$burn_iter + i],
        1,
        unique
      ),
      length
    ) > 20
  )
  if (pmwg_adapted) {
    pmwg_args$adapted = TRUE
    pmwg_args$adapt_iter = pwmg_args$burn_iter + i
    break
  }
}
close(pb)

# Sampling Phase
num_particles <- pmwg_args$sample_particles
cat("Phase 3: Sampling\n")
# create progress bar
pb <- txtProgressBar(
  min = 0,
  max = pmwg_args$burn_iter / progress_update,
  style = 3
)

for (i in 1:pmwg_args$sample_iter) {
  if (i %% progress_update == 0) {
    setTxtProgressBar(pb, i %/% progress_update)
  }

  single_iter <- gen_sample_pars(init, pmwg_args, prior, particles)

  # Sample new particles for random effects.
  if (cpus > 1) {
    tmp <- sfLapply( # nolint
      x = 1:init$S,
      fun = new_sample,
      data = data,
      num_particles = num_particles,
      mu = single_iter$ptm,
      sig2 = single_iter$pts2,
      particles = particles,
      mix_ratio = c(0.1, 0.2, 0.7)
    )
  } else {
    tmp <- lapply(
      X = 1:init$S,
      FUN = new_sample,
      data = data,
      num_particles = num_particles,
      mu = single_iter$ptm,
      sig2 = single_iter$pts2,
      particles = particles,
      mix_ratio = c(0.5, 0.5, 0.0)
    )
  }
  particles <- array(unlist(tmp), dim = dim(particles))

  # Store results.
  init$latent_theta_mu[, , pmwg_args$adapt_iter + i] <- particles # nolint
  init$param_theta_sigma2[, , pmwg_args$adapt_iter + i] <- pts2 # nolint
  init$param_theta_mu[, pmwg_args$adapt_iter + i] <- ptm
}
close(pb)

if (cpus > 1) sfStop() # nolint

save(
  file = restart_file,
  list = c("pts2_inv", "particles", "ptm", "pts2")
)
save.image("data/output/PMwG.RData")
