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
#' @param group_mean An array of starting values for the group means
#' @param group_var An array of starting values for the group variance
#' @param subject_mean An array of starting values for the subject means.
#' @param display_progress Display a progress bar during sampling
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The sampler object but with initial values set for latent_theta_mu
#' @examples
#' sampler <- pmwgs(forstmann, c("b1", "b2", "b3", "A", "v1", "v2", "t0"))
#' sampler <- init(sampler, group_mean=rnorm(7), group_var=diag(rep(0.01, 7)),
#'                 subject_mean=matrix(rnorm(7*19), ncol=19))
#' @export
init.pmwgs <- function(x, group_mean=NULL, group_var=NULL,
                       subject_mean=NULL, display_progress=TRUE, ...) {
  # If no starting point for group mean just use zeros
  if (is.null(group_mean)) group_mean <- stats::rnorm(x$n_pars, sd = 5)
  # If no starting point for group var just sample from inverse wishart
  if (is.null(group_var)) group_var <- MCMCpack::riwish(20, diag(x$n_pars))
  if (is.null(subject_mean)) {
    n_particles <- 1000  #GC: Fixed val here
    subject_mean <- array(NA, dim = c(x$n_pars, x$n_subjects))
    if (display_progress) {
      cat("Sampling Initial values for random effects\n")
      pb <- utils::txtProgressBar(min = 0, max = x$n_subjects, style = 3)
    }
    for (s in 1:x$n_subjects) {
      if (display_progress) utils::setTxtProgressBar(pb, s)
      particles <- mvtnorm::rmvnorm(n_particles, group_mean, group_var)
      colnames(particles) <- rownames(x$samples$group_mean) # preserve par names
      lw <- apply(
        particles,
        1,
        x$llfunc,
        data = x$data[x$data$subject == s, ]
      )
      weight <- exp(lw - max(lw))
      subject_mean[, s] <- particles[
        sample(x = n_particles, size = 1, prob = weight),
      ]
    }
    if (display_progress) close(pb)
  }
  x$init <- TRUE
  x$samples$group_mean[, 1] <- group_mean
  x$samples$group_var[, , 1] <- group_var
  x$samples$subject_mean[, , 1] <- subject_mean
  x
}


#' Initialise variables needed for individual loops within PMwG
#'
#' Takes the pmwgs object and sets up sampling loop variables
#'
#' @param init The list of parameter names to be used in the model
#' @param pmwg_args A list containing arguments to the PMwG sampling, including
#'   number of iterations for each stage.
#' @param prior A list containing priors for variance, means etc.
#' @param particles A multidimensional array containing accepted particles from
#'   the sampling.
#'
#' @return A list of generated variables that can be modified after the fact
#' @examples
#' sampler <- pmwgs(forstmann, c("b1", "b2", "b3", "A", "v1", "v2", "t0"))
#' @export
gen_sample_pars <- function(init, pmwg_args, prior, particles) {
  # Sample population-level parameters.
  var_mu <- MASS::ginv(init$S * pts2_inv + prior$mu_sigma2_inv)
  mean_mu <- as.vector(
    var_mu %*% (pts2_inv %*% apply(particles, 1, sum))
  )
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  # New sample for mu.
  ptm <- mvtnorm::rmvnorm(
    1,
    mean_mu,
    chol_var_mu %*% t(chol_var_mu)
  )[1, ]
  names(ptm) <- init$par_names

  theta_temp <- particles - ptm
  cov_temp <- (theta_temp) %*% (t(theta_temp))
  B_half <- 2 * init$v_half * diag(1 / init$a_half) + cov_temp
  pts2 <- MCMCpack::riwish(init$k_half, B_half) # New sample for sigma.
  pts2_inv <- MASS::ginv(pts2)

  # Sample new mixing weights.
  init$a_half <- 1 / stats::rgamma(
    n = init$num_par,
    shape = init$v_shape,
    scale = 1 / (init$v_half + diag(pts2_inv) + init$A_half)
  )
  list(ptm = ptm, pts2 = pts2, pts2_inv = pts2_inv)
}
