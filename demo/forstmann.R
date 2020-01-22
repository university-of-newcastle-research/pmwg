rm(list = ls())
require(psamplers)
require(rtdists)

# Specify the log likelihood function -----------------------------------------
# This is the likelihood function that we will use for the forstmann dataset.
# It follows the psamplers requirements of having the first element (x) be
# the values for the parameters for which we want to calculate the likelihood.
# It also as it's second argument takes the data.frame object (which is a
# subset of all data just for this subject).
lba_loglike <- function(x, data, sample = FALSE) {
  x <- exp(x)
  if (any(data$rt < x["t0"])) {
    return(-1e10)
  }
  # This is faster than "paste".
  bs <- x["A"] + x[c("b1", "b2", "b3")][data$condition]
  
  if (sample) {
    out <- rtdists::rLBA(n = nrow(data), # nolint
                         A = x["A"],
                         b = bs,
                         t0 = x["t0"],
                         mean_v = x[c("v1", "v2")],
                         sd_v = c(1, 1),
                         distribution = "norm",
                         silent = TRUE)
  } else {
    out <- rtdists::dLBA(rt = data$rt, # nolint
                         response = data$correct,
                         A = x["A"],
                         b = bs,
                         t0 = x["t0"],
                         mean_v = x[c("v1", "v2")],
                         sd_v = c(1, 1),
                         distribution = "norm",
                         silent = TRUE)
    bad <- (out < 1e-10) | (!is.finite(out))
    out[bad] <- 1e-10
    out <- sum(log(out))
  }
  out
}

# Specify the parameters and priors -------------------------------------------

# Vars used for controlling the run
pars <- c("b1", "b2", "b3", "A", "v1", "v2", "t0")
priors <- list(
  theta_mu = rep(0, length(pars)),
  theta_sig = diag(rep(1, length(pars)))
)

# Create the Particle Metropolis within Gibbs sampler object ------------------

sampler <- pmwgs(
  data = forstmann,
  parameters = pars,
  prior = priors
)

start_points <- list(
  mu = c(.2, .2, .2, .4, .3, 1.3, -2),
  sig2 = diag(rep(.01, length(pars)))
)

sampler <- init(sampler, theta_mu = start_points$mu,
                theta_sig = start_points$sig2)

burned <- run_stage(sampler, stage = "burn")

adapted <- run_stage(burned, stage = "adapt")

sampled <- run_stage(adapted, stage = "sample")

save.image("data/output/PMwG.RData")
