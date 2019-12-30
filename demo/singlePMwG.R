library(mvtnorm) ## for the multivariate normal.
library(MASS) ## for matrix inverse.
library(MCMCpack) ## for the inverse wishart random numbers.
library(psamplers)
rm(list = ls())

# Vars used for controlling the run
pars <- c("b1", "b2", "b3", "A", "v1", "v2", "t0")
priors <- list(
  group_mean = rep(0, length(pars)),
  group_var = diag(rep(1, length(pars)))
)

sampler <- pmwgs(
  data = forstmann,
  parameters = pars,
  prior = priors
)

start_points <- list(
  mu = c(.2, .2, .2, .4, .3, 1.3, -2),
  sig2 = diag(rep(.01, length(pars)))
)

sampler <- init(sampler, group_mean = start_points$mu,
                group_var = start_points$sig2)

burned <- run_stage(sampler, stage = "burn")

adapted <- run_stage(burned, stage = "adapt")

sampled <- run_stage(adapted, stage = "sample")

# # Create conditional means/variances
# for (s in 1:init$S) {
#   cparms <- conditional_parms(
#     init,
#     ptm,
#     pts2,
#     s,
#     pmwg_args$burn_iter + 1,
#     pmwg_args$burn_iter + pmwg_args$adapt_iter
#   )
#   proposal_means[, s] <- cparms$cmeans
#   proposal_sigmas[, , s] <- cparms$cvars
# }
# 
# # Sampling Phase
# num_particles <- pmwg_args$sample_particles
# cat("Phase 3: Sampling\n")
# # create progress bar
# pb <- txtProgressBar(
#   min = 0,
#   max = pmwg_args$sample_iter / progress_update,
#   style = 3
# )
# 
# for (i in 1:pmwg_args$sample_iter) {
#   if (i %% progress_update == 0) {
#     setTxtProgressBar(pb, i %/% progress_update)
#   }
# 
#   single_iter <- gen_sample_pars(init, pmwg_args, prior, particles)
#   ptm <- single_iter$ptm
#   pts2 <- single_iter$pts2
# 
#   # Sample new particles for random effects.
#   tmp <- lapply(
#     X = 1:init$S,
#     FUN = new_sample,
#     data = data,
#     num_particles = num_particles,
#     mu = ptm,
#     sig2 = pts2,
#     particles = particles,
#     mix_ratio = c(0.1, 0.2, 0.7)
#   )
#   particles <- array(unlist(tmp), dim = dim(particles))
# 
#   # Store results.
#   init$latent_theta_mu[, , pmwg_args$adapt_iter + i] <- particles # nolint
#   init$param_theta_sigma2[, , pmwg_args$adapt_iter + i] <- pts2 # nolint
#   init$param_theta_mu[, pmwg_args$adapt_iter + i] <- ptm
# }
# close(pb)

save.image("data/output/PMwG.RData")
