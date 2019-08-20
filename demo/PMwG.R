library(rtdists)
library(mvtnorm) ## for the multivariate normal.
library(MASS) ## for matrix inverse.
library(MCMCpack) ## for the inverse wishart random numbers.
library(psamplers)
rm(list = ls())

runner_args <- list("restart" = FALSE,
                    "cpus" = 1,
                    "restart_file" = "data/output/restart.RData")
n.particles <- 100
n.iterations <- 100


# Load Forstmann et al.'s data.
data <- read.csv("data/data.csv", header = FALSE)
names(data) <- c("subject", "rt", "correct", "condition")
S <- 3

parameter.names <- c("b1", "b2", "b3", "A", "v1", "v2", "t0")
n.parameters <- length(parameter.names)
start_points_mu <- c(.2, .2, .2, .4, .3, 1.3, -2)
start_points_sig2 <- diag(rep(.01, n.parameters))

# Tuning settings for the Gibbs steps.
v.half <- 2
A.half <- 1

# Storage for the samples.
latent.theta.mu <- array(NA, dim = c(n.parameters, S, n.iterations), dimnames = list(parameter.names, NULL, NULL))
param.theta.mu <- latent.theta.mu[, 1, ]
param.theta.sigma2 <- array(NA, dim = c(n.parameters, n.parameters, n.iterations), dimnames = list(parameter.names, parameter.names, NULL))

# Make single-iteration-sized versions, for easier reading of code below.
ptm <- param.theta.mu[, 1]
pts2 <- param.theta.sigma2[, , 1]

# Start points for the population-level parameters only. Hard coded here, just
# for convenience.
ptm[1:n.parameters] <- start_points_mu # Weird subscripts maintains the naming.
pts2 <- start_points_sig2  # Who knows??
pts2.inv <- ginv(pts2) # Because this is calculated near the end of the main loop, needs initialising for iter=1.

# Priors.
prior.mu.mean <- rep(0, n.parameters)
prior.mu.sigma2 <- diag(rep(1, n.parameters))


# Things I save rather than re-compute inside the loops.
k.half <- v.half + n.parameters - 1 + S
v.shape <- (v.half + n.parameters) / 2
prior.mu.sigma2.inv <- ginv(prior.mu.sigma2)

# Sample the initial values for the random effects. Algorithm is same
# as for the main resampling down below.
particles <- array(dim = c(length(ptm), S))
for (s in 1:S) {
  proposals <- rmvnorm(n.particles, ptm, pts2)
  lw <- apply(proposals, 1, lba_loglike, data = data[data$subject == s, ])
  weight <- exp(lw - max(lw))
  particles[, s] <- proposals[sample(x = n.particles, size = 1, prob = weight), ]
}

# Sample the mixture variables' initial values.
a.half <- 1 / rgamma(n = n.parameters, shape = 0.5, scale = 1)

if (runner_args$restart) {
  cat("\nRestarting from saved run.\n")
  load(runner_args$restart_file)
  # Check a couple of things.
  if (dim(particles)[1] != n.parameters) stop("Restart does not match size (params).")
  if (dim(particles)[2] != S) stop("Restart does not match size (subjects).")
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
    "param.theta.mu",
    "param.theta.sig2",
    "latent.theta.mu"
  ))
}

total <- 20
# create progress bar
pb <- txtProgressBar(min = 0, max = total, style = 3)
for(i in 1:total){
   Sys.sleep(0.1)
   # update progress bar
   setTxtProgressBar(pb, i)
}
close(pb)

for (i in 1:n.iterations) {
  cat("\t", i)
  # Sample population-level parameters.
  var_mu <- ginv(S * pts2.inv + prior.mu.sigma2.inv)
  mean_mu <- as.vector(var_mu %*% (pts2.inv %*% apply(particles, 1, sum)))
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  ptm <- rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ] # New sample for mu.
  names(ptm) <- parameter.names

  theta.temp <- particles - ptm
  # cov.temp=array(0,dim=c(n.parameters,n.parameters))
  # for (j in 1:S) cov.temp=cov.temp+(theta.temp[,j])%*%t(theta.temp[,j])
  cov.temp <- (theta.temp) %*% (t(theta.temp))
  B.half <- 2 * v.half * diag(1 / a.half) + cov.temp
  pts2 <- riwish(k.half, B.half) # New sample for sigma.
  pts2.inv <- ginv(pts2)

  # Sample new mixing weights.
  a.half <- 1 / rgamma(n = n.parameters, shape = v.shape, scale = 1 / (v.half + diag(pts2.inv) + A.half))

  # Sample new particles for random effects.
  if (runner_args$cpus > 1) {
    tmp <- sfLapply(x = 1:S, fun = new_sample, data = data, num_proposals = n.particles, mu = ptm, sig2 = pts2, particles = particles)
  } else {
    tmp <- lapply(X = 1:S, FUN = new_sample, data = data, num_proposals = n.particles, mu = ptm, sig2 = pts2, particles = particles)
  }
  particles <- array(unlist(tmp), dim = dim(particles))

  # Store results.
  latent.theta.mu[, , i] <- particles
  param.theta.sigma2[, , i] <- pts2
  param.theta.mu[, i] <- ptm
}
if (runner_args$cpus > 1) sfStop()

save(file = runner_args$restart_file, list = c("pts2.inv", "particles", "ptm", "pts2"))
save.image("data/output/PMwG.RData")
