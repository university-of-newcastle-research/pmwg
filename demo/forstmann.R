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

save.image("data/output/PMwG.RData")
