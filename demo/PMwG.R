library(rtdists)
library(mvtnorm) ## for the multivariate normal.
library(MASS) ## for matrix inverse.
library(MCMCpack) ## for the inverse wishart random numbers.
library(psamplers)
rm(list = ls())

runner_args <- list(
  "restart" = FALSE,
  "cpus" = 1,
  "restart_file" = "data/output/restart.RData",
  "progress_update" = 10
)

pmwg_args <- list(
  "adaptation_particles" = 100,
  "burnin_iterations" = 500,
  "sampling_particles" = 100,
  "sampling_iterations" = 1000,
  "max_iterations" = 10000,
  "thin" = 1
)
# pmwg_args <- list(
#   "adaptation_particles" = 1000,
#   "burnin_iterations" = 500,
#   "sampling_particles" = 100,
#   "sampling_iterations" = 10000,
#   "max_iterations" = 100000,
#   "thin" = 1
# )

# Load Forstmann et al.'s data.
data <- read.csv("data/data.csv", header = FALSE)
names(data) <- c("subject", "rt", "correct", "condition")
S <- length(unique(data$subject))

parameters <- c("b1", "b2", "b3", "A", "v1", "v2", "t0")
num_parameters <- length(parameters)
start_points_mu <- c(.2, .2, .2, .4, .3, 1.3, -2)
start_points_sig2 <- diag(rep(.01, num_parameters))

# Tuning settings for the Gibbs steps
v_half <- 2
A_half <- 1

# theta is the parameter values, mu is mean of normnal distribution and sigma2 is variance
# Storage for the samples.
latent_theta_mu <- array(
  NA,
  dim = c(num_parameters, S, pmwg_args$sampling_iterations),
  dimnames = list(parameters, NULL, NULL)
)
param_theta_mu <- latent_theta_mu[, 1, ]
param_theta_sigma2 <- array(
  NA,
  dim = c(num_parameters, num_parameters, pmwg_args$sampling_iterations),
  dimnames = list(parameters, parameters, NULL)
)

# Make single-iteration-sized versions, for easier reading of code below.
ptm <- param_theta_mu[, 1]
pts2 <- param_theta_sigma2[, , 1] # nolint

# Start points for the population-level parameters only. Hard coded here, just
# for convenience.
ptm[1:num_parameters] <- start_points_mu # Weird subscripts maintains the naming.
pts2 <- start_points_sig2 # Who knows??
# Because this is calculated near the end of the main loop, needs initialising for iter=1.
pts2_inv <- ginv(pts2)

# Priors.
prior_mu_mean <- rep(0, num_parameters)
prior_mu_sigma2 <- diag(rep(1, num_parameters))


# Things I save rather than re-compute inside the loops.
k_half <- v_half + num_parameters - 1 + S
v_shape <- (v_half + num_parameters) / 2
prior_mu_sigma2_inv <- ginv(prior_mu_sigma2)

# Sample the initial values for the random effects. Algorithm is same
# as for the main resampling down below.
particles <- array(dim = c(length(ptm), S))
num_particles <- pmwg_args$adaptation_particles
for (s in 1:S) {
  proposals <- rmvnorm(num_particles, ptm, pts2)
  lw <- apply(proposals, 1, lba_loglike, data = data[data$subject == s, ])
  weight <- exp(lw - max(lw))
  particles[, s] <- proposals[
    sample(x = num_particles, size = 1, prob = weight),
  ]
}

# Sample the mixture variables' initial values.
a_half <- 1 / rgamma(n = num_parameters, shape = 0.5, scale = 1)

if (runner_args$restart) {
  cat("\nRestarting from saved run.\n")
  load(runner_args$restart_file)
  # Check a couple of things.
  if ((dim(particles)[1] != num_parameters) || (dim(particles)[2] != S)) {
    stop("Restart does not match size (subjects).")
  }
}

if (runner_args$cpus > 1) {
  # You need the suggested package for this function
  library(snowfall)

  sfInit(parallel = TRUE, cpus = runner_args$cpus)
  sfClusterSetupRNG()
  sfLibrary(rtdists)
  sfLibrary(mvtnorm)
  sfLibrary(MASS) # For matrix inverse.
  sfLibrary(MCMCpack) # For the inverse Wishart random numbers.

  sfExportAll(except = c(
    "param_theta_mu",
    "param_theta_sig2",
    "latent_theta_mu"
  ))
}

# create progress bar
pu <- runner_args$progress_update
pb <- txtProgressBar(min = 0, max = pmwg_args$sampling_iterations/pu, style = 3)

for (i in 1:pmwg_args$sampling_iterations) {
  if (i %% pu == 0){
    setTxtProgressBar(pb, i %/% pu)
  }
  # Sample population-level parameters.
  var_mu <- ginv(S * pts2_inv + prior_mu_sigma2_inv)
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

save(file = runner_args$restart_file, list = c("pts2_inv", "particles", "ptm", "pts2"))
save.image("data/output/PMwG.RData")
