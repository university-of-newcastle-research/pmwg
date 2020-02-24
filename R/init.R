#' Initialise values for the random effects
#'
#' Takes the starting values for the group mean and variance, and subject level
#' means. All arrays must match the appropriate shape.
#'
#' For example, with 5 parameters and 10 subjects, the group means must be a
#' vector of length 5, the group variance must be an array of 5 x 5, and the
#' subject means must be 5 x 10.
#'
#' Alternatively the if argument values for the starting points are left at the
#' default (NULL) then starting points will be sampled from the prior for group
#' level values, or for subject level means they will be sampled from the
#' multivariate normal using the group level means and variance.
#'
#' @param x The sampler object that provides the parameters.
#' @param theta_mu An array of starting values for the group means
#' @param theta_sig An array of starting values for the group covariance matrix
#' @param display_progress Display a progress bar during sampling
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The sampler object but with initial values set for latent_theta_mu
#' @examples
#' lba_ll <- function(x, data) {
#'   sum(
#'     log(
#'       rtdists::dLBA(rt = data$rt,
#'                     response = data$correct,
#'                     A = x["A"],
#'                     b = bs,
#'                     t0 = x["t0"],
#'                     mean_v = x[c("v1", "v2")],
#'                     sd_v = c(1, 1),
#'                     silent = TRUE)
#'     )
#'   )
#' }
#' sampler <- pmwgs(forstmann,
#'                  c("b1", "b2", "b3", "A", "v1", "v2", "t0"),
#'                  lba_ll
#'            )
#' sampler <- init(sampler, theta_mu=rnorm(7), theta_sig=diag(rep(0.01, 7)))
#' @export
init.pmwgs <- function(x, theta_mu=NULL, theta_sig=NULL,
                       display_progress=TRUE, ...) {
  # If no starting point for group mean just use zeros
  if (is.null(theta_mu)) theta_mu <- stats::rnorm(x$n_pars, sd = 1)
  # If no starting point for group var just sample from inverse wishart
  if (is.null(theta_sig)) theta_sig <- MCMCpack::riwish(20, diag(x$n_pars))
  n_particles <- 1000  #GC: Fixed val here
  alpha <- array(NA, dim = c(x$n_pars, x$n_subjects))
  if (display_progress) {
    cat("Sampling Initial values for random effects\n")
    pb <- utils::txtProgressBar(min = 0, max = x$n_subjects, style = 3)
  }
  likelihoods <- array(NA_real_, dim = c(x$n_subjects))
  for (s in 1:x$n_subjects) {
    if (display_progress) utils::setTxtProgressBar(pb, s)
    particles <- mvtnorm::rmvnorm(n_particles, theta_mu, theta_sig)
    colnames(particles) <- rownames(x$samples$theta_mu) # preserve par names
    lw <- apply(
      particles,
      1,
      x$ll_func,
      data = x$data[x$data$subject == x$subjects[s], ]
    )
    weight <- exp(lw - max(lw))
    idx <- sample(x = n_particles, size = 1, prob = weight)
    alpha[, s] <- particles[idx, ]
    likelihoods[s] <- lw[idx]
  }
  if (display_progress) close(pb)
  x$init <- TRUE
  x$samples$theta_mu[, 1] <- theta_mu
  x$samples$theta_sig[, , 1] <- theta_sig
  x$samples$alpha[, , 1] <- alpha
  x$samples$last_theta_sig_inverse <- MASS::ginv(theta_sig)
  x$samples$subj_ll[, 1] <- likelihoods
  x$samples$idx <- 1
  x
}


#' Initialise variables needed for individual loops within PMwG
#'
#' Takes the pmwgs object and sets up sampling loop variables
#'
#' @param samples The list containing the samples from the current run, or from
#'   the master storage in the sampler
#' @param sampler The pmwgs object from which to generate the new group
#'   parameters.
#'
#' @return A list of generated variables that can be modified after the fact
#' @examples
#' # No example yet
#' @export
new_group_pars <- function(samples, sampler) {
  # Get single iter versions, gm = theta_mu, gv = theta_sig
  last <- last_sample(samples)
  hyper <- attributes(sampler)

  # Here mu is group mean, so we are getting mean and variance
  var_mu <- MASS::ginv(
    sampler$n_subjects * last$gvi + sampler$prior$theta_sig_inv
  )
  mean_mu <- as.vector(var_mu %*% (last$gvi %*% apply(last$sm, 1, sum)))
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  # New sample for mu.
  gm <- mvtnorm::rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
  names(gm) <- sampler$par_names

  # New values for group var
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
  list(gm = gm, gv = gv, gvi = gvi, a_half = a_half, sm = last$sm)
}
