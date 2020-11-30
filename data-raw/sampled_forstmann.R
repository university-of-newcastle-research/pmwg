library(pmwg)
library(rtdists)
set.seed(24061795) # Weber birthday

# Specify the log likelihood function
lba_loglike <- function(x, data) {
  x <- exp(x)
  if (any(data$rt < x["t0"])) {
    return(-1e10)
  }
  # This is faster than "paste".
  bs <- x["A"] + x[c("b1", "b2", "b3")][data$condition]

  out <- rtdists::dLBA(
    rt = data$rt, # nolint
    response = data$correct,
    A = x["A"],
    b = bs,
    t0 = x["t0"],
    mean_v = x[c("v1", "v2")],
    sd_v = c(1, 1),
    distribution = "norm",
    silent = TRUE
  )
  bad <- (out < 1e-10) | (!is.finite(out))
  out[bad] <- 1e-10
  out <- sum(log(out))
  out
}

# Add a correct column for our log likelihood
forstmann$correct <- (forstmann$stim == forstmann$resp) + 1

# Create the sampler object with your data, parameter names, loglike function and
# A list of priors for theta_mu_mean (mean of model parameters) and theta_mu_var (covariance of model parameters)
sampler <- pmwgs(
  data = forstmann,
  pars = c("b1", "b2", "b3", "A", "v1", "v2", "t0"),
  ll_func = lba_loglike,
)

# Initialise (generate first random effects for sampler)
sampler <- init(sampler)  # Can also pass start points for sampler

# Run each stage of the sampler, can adjust number of particles on each
sampler <- run_stage(sampler, iter=20, particles = 100, stage="burn")
sampler <- run_stage(sampler, stage="adapt", particles = 100, iter=100)
sampler <- run_stage(sampler, particles=50, iter=80, stage="sample")

sampled_forstmann <- sampler
sampled_forstmann$data <- NULL

usethis::use_data(sampled_forstmann, overwrite = TRUE)
