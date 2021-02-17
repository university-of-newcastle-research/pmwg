devtools::load_all()
library(rtdists)
set.seed(24061795) # Weber birthday

# Specify the log likelihood function
fast_lba_ll3b <- function(x, data, sample = FALSE) {
  x <- exp(x)
  if (any(data$rt < x["t0"])) {
    return(-1e10)
  }

  out <- numeric(nrow(data))

  all_b <- numeric(nrow(data))
  vlist <- list(
    "v.1" = numeric(nrow(data)),
    "v.2" = numeric(nrow(data))
  )
  stim <- levels(data$stim)
  con <- levels(data$condition)

  for (c in con) {
    for (s in stim) {
      use <- data$condition == c & data$stim == s
      if (any(use)) {
        bs <- x[paste0("b.", c)] + x["A"]
        all_b[use] <- bs
        vc <- x["vc"]
        ve <- x["ve"]
        if (s == 1) {
          vlist$v.1[use] <- vc
          vlist$v.2[use] <- ve
        } else {
          vlist$v.1[use] <- ve
          vlist$v.2[use] <- vc
        }
      }
    }
  }

  out <- dLBA(
    rt = data$rt,
    response = data$resp,
    A = x["A"],
    b = all_b,
    mean_v = vlist,
    sd_v = c(1, 1),
    t0 = x["t0"],
    distribution = "norm",
    silent = TRUE
  )

  bad <- (out < 1e-10) | (!is.finite(out))
  out[bad] <- 1e-10
  out <- sum(log(out))
  return(out)
}

# Create the sampler object with your data, parameter names, loglike function
# and a list of priors for theta_mu_mean (mean of model parameters) and
# theta_mu_var (covariance of model parameters)
sampler <- pmwgs(
  data = forstmann,
  pars <- c("b1", "b2", "b3", "A", "ve", "vc", "t0"),
  ll_func = fast_lba_ll3b,
)

# Initialise (generate first random effects for sampler)
sampler <- init(sampler) # Can also pass start points for sampler

# Run each stage of the sampler, can adjust number of particles on each
sampler <- run_stage(sampler, iter = 20, particles = 100, stage = "burn")
sampler <- run_stage(sampler, stage = "adapt", particles = 100, iter = 100)
sampler <- run_stage(sampler, particles = 50, iter = 80, stage = "sample")

sampled_forstmann <- sampler
# Strip out the data element for more efficient storage
sampled_forstmann$data <- NULL

usethis::use_data(sampled_forstmann, overwrite = TRUE)
