# Load libs, data ---------------------------------------------------------
# Specify log likelihood func (wrapper on rtdists)-------------------------
# Specify the parameters and priors ---------------------------------------
pars <- c("b1", "b2", "b3", "A", "v1", "v2", "t0")
priors <- list(
  theta_mu = rep(0, length(pars)),
  theta_sig = diag(rep(1, length(pars)))
)
# Create the Particle Metropolis within Gibbs sampler object --------------
sampler <- pmwgs(
  data = forstmann,
  pars = pars,
  ll_func = lba_loglike,
  prior = priors
)
# Start points ------------------------------------------------------------
start_points <- list(
  mu = c(.2, .2, .2, .4, .3, 1.3, -2),
  sig2 = diag(rep(.01, length(pars)))
)
# No setup yet for number iterations etc - all in run_stage ---------------







# Initialise the sampler --------------------------------------------------
sampler <- init(sampler, theta_mu = start_points$mu,
                theta_sig = start_points$sig2)
# Copy of init code -------------------------------------------------------
init.pmwgs <- function(x, theta_mu=NULL, theta_sig=NULL, display_progress=TRUE, ...) {
  # Setup num particles, creating output array alpha etc
  for (s in 1:x$n_subjects) {
    particles <- mvtnorm::rmvnorm(n_particles, theta_mu, theta_sig)
    #Using rtdists so not needed to manipulate particles into diff format
    lw <- apply(particles, 1, x$ll_func, data = x$data[x$data$subject == x$subjects[s], ])
    # Our loglike function tests for impossible numbers etc not needed here
    weight <- exp(lw - max(lw))

    # We don't divide by sum because sample function can take weights not summing to 1
    # Additional check for negative weights, etc
    idx <- sample(x = n_particles, size = 1, prob = weight)
    alpha[, s] <- particles[idx, ]
    likelihoods[s] <- lw[idx]
  }
  # Set stores in x (pmwgs object) with alpha, likelihoods etc
}
# No more setup needed - ready to run stages-------------------------------
# Run_stage function, meat of processing ----------------------------------
burned <- run_stage(sampler, stage = "burn")
run_stage.pmwgs <- function(x, stage, iter = 1000, particles = 1000, display_progress = TRUE, n_cores = 1, ...) {
# Select stage, check args build sample storage----------------------------
# If stage is sample, try to create_efficient proposal args ---------------
  for (i in 1:iter) {
# Call Function for gibbs step, plus function contents --------------------
      pars <- new_group_pars(store, x),
new_group_pars <- function(samples, sampler) {
# Get last sample, hyper pars ---------------------------------------------
  var_mu <- MASS::ginv(sampler$n_subjects * last$gvi + sampler$prior$theta_sig_inv)
  mean_mu <- as.vector(var_mu %*% (last$gvi %*% apply(last$sm, 1, sum)))
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  gm <- mvtnorm::rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
  # k_half calculated earlier
  theta_temp <- last$sm - gm
  cov_temp <- (theta_temp) %*% (t(theta_temp))


  
  B_half <- 2 * hyper$v_half * diag(1 / hyper$a_half) + cov_temp #nolint
  gv <- MCMCpack::riwish(hyper$k_half, B_half) # New sample for group variance
  gvi <- MASS::ginv(gv)

  # Sample new mixing weights.
  a_half <- 1 / stats::rgamma(
    n = sampler$n_pars,
    shape = hyper$v_shape,
    scale = 1 / (hyper$v_half + diag(gvi) + hyper$A_half)
  )
#  return calculated values
}
# Set arguments and functions for getting new samples---------------------
# Contents of create_efficient (called earlier) repro'ed for comparison ----
create_efficient <- function(x) {
  # Create storage for conditional pars, by subject run conditional_parms
}
conditional_parms <- function(s, samples) {
  all_samples <- rbind(
    samples$alpha[, s, ],
    samples$theta_mu[, ],
    pts2_unwound
  ) # pts_unwound is:
  #y <- t(chol(var_matrix))  # take cholesky factor (top right), transpose
  #diag(y) <- log(diag(y))   # log the diagonal
  #y[lower.tri(y, diag = TRUE)]  
# Create sample mean sigma for joint random effect and pars ---------------------
  mu_tilde <- apply(all_samples, 1, mean)
  sigma_tilde <- stats::var(t(all_samples))
  #var will computer covariance matrix, we don't correct for non pos definite
# R code calculates conidtional mean/var here, code later for comparison)-----
  condmvn <- condMVNorm::condMVN(...)
# Switching num particles done at higher level----------------------------------------------
# Obtain a new sample -----------------------------------------------------
    tmp <- do.call(apply_fn, fn_args)
new_sample <- function(s, data, num_particles, parameters,
                       efficient_mu = NULL, efficient_sig2 = NULL,
                       mix_ratio = c(0.5, 0.5, 0.0),
                       likelihood_func = NULL,
                       epsilon = 1, subjects = NULL) {
# Already at subject level
# Check args, extract efficient proposal means/last sample pars -----------
# no equivalent here but copied from above---------------------------------
  condmvn <- condMVNorm::condMVN(
    mean = mu_tilde,
    sigma = sigma_tilde,
    dependent.ind = 1:n_par,
    given.ind = (n_par + 1):length(mu_tilde),
    X.given = c(samples$theta_mu[, n_iter],
                unwind(samples$theta_sig[, , n_iter]))
  )
  list(cmeans = condmvn$condMean, cvars = condmvn$condVar)
}
# Create proposals for new particles --------------------------------------
  proposals <- gen_particles(num_particles, mu, sig2, subj_mu, mix_ratio = mix_ratio, prop_mu = e_mu, prop_sig2 = e_sig2, epsilon = epsilon )
  # Except here:
  particle_numbers <- numbers_from_ratio(mix_ratio, num_particles)
  pop_particles <- particle_draws(particle_numbers[1], mu, sig2)
  ind_particles <- particle_draws(particle_numbers[2], particle, sig2 * epsilon)
  eff_particles <- particle_draws(particle_numbers[3], prop_mu, prop_sig2)
  particles <- rbind(pop_particles, ind_particles, eff_particles)
  # Put the current particle in slot 1.
  proposals[1, ] <- subj_mu
# No rearranging necessary -----------------------------------------------

# Again no rearrangin necessary ------------------------------------------
# computing the log density of the LBA given the particles of random effects-----
  lw <- apply(proposals, 1, likelihood_func, data = data[data$subject == subjects[s], ])


  
# Density of random effects proposal given population-level distribution.----
  lp <- mvtnorm::dmvnorm(x = proposals, mean = mu, sigma = sig2, log = TRUE)
  # Density of proposals given proposal distribution.
  prop_density <- mvtnorm::dmvnorm(x = proposals, mean = subj_mu, sigma = sig2)
  # Density of efficient proposals
  if (mix_ratio[3] != 0) {
    eff_density <- mvtnorm::dmvnorm(
      x = proposals,
      mean = e_mu,
      sigma = e_sig2
    )
  } else {
    eff_density <- 0
  }
l
  lm <- log(mix_ratio[1] * exp(lp) +
    (mix_ratio[2] * prop_density) +
    (mix_ratio[3] * eff_density))
  # log of importance weights.
  l <- lw + lp - lm
# No checks here for imaginary number of logw ----------------------------
  weights <- exp(l - max(l))





# Sample winning particle ------------------------------------------------
  idx <- sample(x = num_particles, size = 1, prob = weights)
  winner <- proposals[idx, ]
  attr(winner, "ll") <- lw[idx]
  winner
}





# storing the MCMC draws------------------------------------------------
# Check for number unique values in random effects----------------------
# Save up to the user (RData for example)-------------------------------
