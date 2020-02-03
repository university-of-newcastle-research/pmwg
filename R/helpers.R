#' Unwinds variance matrix to a vector
#'
#' Takes a variance matrix and unwind to a vector via Cholesky then log
#'
#' @param var_matrix A variance matrix
#'
#' @return The unwound matrix as a vector
#' @examples
#' psamplers:::unwind(diag(rep(1, 7)))
#' @keywords internal
unwind <- function(var_matrix, ...) {
  y <- t(chol(var_matrix))
  diag(y) <- log(diag(y))
  y[lower.tri(y, diag = TRUE)]
}

#' Winds a variance vector back to a vector
#'
#' The reverse of the unwind function, takes a variance vector and windows back into matrix
#'
#' @param var_vector A variance vector
#'
#' @return The wound vector as a matrix
#' @examples
#' psamplers:::wind(diag(rep(1, 7)))
#' @keywords internal
wind <- function(var_vector, ...) {
  n <- sqrt(2 * length(var_vector) + 0.25) - 0.5 ## Dim of matrix.
  if ((n * n + n) != (2 * length(var_vector))) stop("Wrong sizes in unwind.")
  out <- array(0, dim = c(n, n))
  out[lower.tri(out, diag = TRUE)] <- var_vector
  diag(out) <- exp(diag(out))
  out %*% t(out)
}


#' Check and normalise the number of each particle type from the mix_ratio
#'
#' Takes a mix ratio vector (3 x float) and a number of particles to generate and returns
#' a vector containing the number of each particle type to generate
#'
#' @param mix_ratio A vector of floats betwen 0 and 1 and summing to 1 which give the ratio
#'   of particles to generate from the population level parameters, the individual random
#'   effects and the conditional parameters repectively
#' @param num_particles The total number of particles to generate using a combination of the
#'   three methods.
#'
#' @return The wound vector as a matrix
#' @examples
#' psamplers:::numbers_from_ratio(c(0.1, 0.3, 0.6))
#' @keywords internal
numbers_from_ratio <- function(mix_ratio, num_particles = 1000) {
  if (!isTRUE(all.equal(sum(mix_ratio), 1))) {
    stop("The elements of the mix_ratio vector must sum to 1")
  }
  if (length(mix_ratio) != 3) {
    stop("mix_ratio vector must have three elements which sum to 1")
  }
  numbers <- stats::rbinom(n = 2, size = num_particles, prob = mix_ratio)
  if (mix_ratio[3] == 0) {
    numbers[3] <- 0
    numbers[2] <- num_particles - numbers[1]
  } else {
    numbers[3] <- num_particles - sum(numbers)
  }
  numbers
}


#' Check for efficient proposals if necessary
#'
#' Takes a mix ratio vector (3 x float) and the efficient proposal mu and sigma
#' If efficient proposals are to be used (mix_ratio[3] > 0) then test the
#' efficient proposal values to see whether they are not null and appropriate.
#'
#' @param efficient_mu The mu value for the efficient proposals
#' @param efficient_sig2 The sigma value for the efficient proposals
#' @param mix_ratio A vector of floats betwen 0 and 1 and summing to 1 which give the ratio
#'   of particles to generate from the population level parameters, the individual random
#'   effects and the conditional parameters repectively
#'
#' @return nothing, stops operation on incorrect combiation of parameters.
#' @examples
#' psamplers:::check_efficient(c(0.1, 0.9, 0.0), NULL, NULL)
#' @keywords internal
check_efficient <- function(mix_ratio, efficient_mu, efficient_sig2) {
  if (mix_ratio[3] != 0) {
    if (is.null(efficient_mu) || is.null(efficient_sig2)) {
      stop(
        paste0(
          "Mu and sigma from efficient conditional ",
          "proposals must be provided for mix_ratio[3] > 0"
        )
      )
    }
  }
}


#' Extract relevant samples from the list for conditional dist calc
#'
#' From the existing samples, extract relevant samples for the creation of
#' the proposal distribution.
#'
#' @param samples The samples list containing all samples from the pmwgs object
#'
#' @return A list containing only appopriate samples (non init/burnin samples)
#' @examples
#' # No example yet
#' @keywords internal
extract_samples <- function(samples) {
  sample_filter <- samples$stage %in% c("adapt", "sample")
  list(
    theta_mu = samples$theta_mu[, sample_filter],
    theta_sig = samples$theta_sig[, , sample_filter],
    alpha = samples$alpha[, , sample_filter]
  )
}


#' Create distribution parameters for efficient proposals
#'
#' From the existing samples, create a proposal distribution for drawing
#' efficient samples from.
#'
#' @param x The current pmwgs object
#'
#' @return A list containing the mu and sigma for the proposal distribution.            
#' @examples
#' # No example yet
#' @keywords internal
create_efficient <- function(x) {
  proposal_means <- array(dim = c(x$n_pars, x$n_subjects))
  proposal_sigmas <- array(dim = c(x$n_pars, x$n_pars, x$n_subjects))
  for (s in 1:x$n_subjects) {
    cparms <- conditional_parms(
      extract_samples(x$samples),
      s
    )
    proposal_means[, s] <- cparms$cmeans
    proposal_sigmas[, , s] <- cparms$cvars
  }
  list(
    efficient_mu = proposal_means,
    efficient_sig2 = proposal_sigmas
  )
}


#' Generate a cloud of particles from a multivariate normal distribution
#'
#' Takes the mean and variance for a multivariate normal distribution, as well as the number of
#' particles to generate and return random draws from the multivariate normal if the numbers of
#' particles is > 0, otherwise return NULL. At least one of mean or sigma must be provided.
#'
#' @param n number of observations
#' @param mu mean vector
#' @param covar covariance matrix
#'
#' @return If n > 0 returns n draws from the multivariate normal with mean and sigma, otherwise returns NULL
#' @examples
#' psamplers:::particle_draws(100, rep(0.2, 7), diag(rep(7)))
#' psamplers:::particle_draws(0, rep(0.2, 7), diag(rep(7)))
#' @keywords internal
particle_draws <- function(n, mu, covar) {
  if (n <= 0) {
    return(NULL)
  }
  mvtnorm::rmvnorm(n, mean = mu, sigma = covar)
}


#' Obtain the efficent mu and sigma from the adaptation phase draws
#'
#' @param samples A list containing previous samples
#' @param s current subject number
#'
#' @return A list containing the conditional mean and variances for this subject
#' @examples
#' # No example yet
#' @export
conditional_parms <- function(samples, s) {
  gmdim <- dim(samples$theta_mu)
  n_par <- gmdim[1]
  n_iter <- gmdim[2]
  pts2_unwound <- apply(
    samples$theta_sig,
    3,
    unwind
  )
  all_samples <- rbind(
    samples$alpha[, s, ],
    samples$theta_mu[, ],
    pts2_unwound
  )
  mu_tilde <- apply(all_samples, 1, mean)
  sigma_tilde <- stats::var(t(all_samples))
  condmvn <- condMVNorm::condMVN(
    mean = mu_tilde,
    sigma = sigma_tilde,
    dependent.ind = 1:n_par,
    given.ind = (n_par + 1):length(mu_tilde),
    # GC: Note, not sure what is happening here:v (Was ptm/pts2 now last sample)
    X.given = c(samples$theta_mu[, n_iter],
                unwind(samples$theta_sig[, , n_iter]))
  )
  list(cmeans = condmvn$condMean, cvars = condmvn$condVar)
}


#' Create a new list for storage samples in the pmwgs object
#'
#' @param par_names The names of each parameter as a character vector
#' @param n_subjects The number of subjects for the subject mean storage.
#' @param iters The number of iterations to be pre-allocated
#' @param stage The stage for which the samples will be created. Should be one
#'   of `c("init", "burn", "adapt", "sample")`
#'
#' @return A list containing the conditional mean and variances for this subject
#' @examples
#' # No example yet
#' @export
sample_store <- function(par_names, n_subjects, iters = 1, stage = "init") {
  n_pars <- length(par_names)
  list(
    alpha = array(
      NA_real_,
      dim = c(n_pars, n_subjects, iters),
      dimnames = list(par_names, NULL, NULL)
    ),
    theta_mu = array(
      NA_real_,
      dim = c(n_pars, iters),
      dimnames = list(par_names, NULL)
    ),
    theta_sig = array(
      NA_real_,
      dim = c(n_pars, n_pars, iters),
      dimnames = list(par_names, par_names, NULL)
    ),
    stage = array(stage, iters)
  )
}


#' Create a list with the last samples in the pmwgs object
#'
#' @param store The list containing samples from t=which to grab the last.
#'
#' @return A list containing the last sample of group mean and variance and
#'   subject means.
#' @examples
#' # No example yet
#' @export
last_sample <- function(store) {
  list(
    gm = store$theta_mu[, store$idx],
    gv = store$theta_sig[, , store$idx],
    sm = store$alpha[, , store$idx],
    gvi = store$last_theta_sig_inv
  )
}


#' Update the main data store with the results of the last stage
#'
#' @param sampler The pmwgs object that we are adding the new samples to
#' @param store The sample storage stage just run
#'
#' @return The pmwgs object with the new samples concatenated to the old
#' @examples
#' # No example yet
#' @export
update_sampler <- function(sampler, store) {
  old_gm <- sampler$samples$theta_mu
  old_gv <- sampler$samples$theta_sig
  old_sm <- sampler$samples$alpha
  old_stage <- sampler$samples$stage
  li <- store$idx

  sampler$samples$theta_mu <- array(c(old_gm, store$theta_mu[, 1:li]),
                                      dim = dim(old_gm) + c(0, li))
  sampler$samples$theta_sig <- array(c(old_gv, store$theta_sig[, , 1:li]),
                                     dim = dim(old_gv) + c(0, 0, li))
  sampler$samples$alpha <- array(c(old_sm, store$alpha[, , 1:li]),
                                        dim = dim(old_sm) + c(0, 0, li))
  sampler$samples$idx <- ncol(sampler$samples$theta_mu)
  sampler$samples$last_theta_sig_inv <- store$last_theta_sig_inv
  sampler$samples$stage <- c(old_stage, store$stage[1:li])
  sampler
}


#' Check whether the adaptation phase has successfully completed
#'
#' @param samples The subject mean samples with which we are working
#' @param unq_vals The number of unique values for each subject
#'
#' @return A boolean TRUE or FALSE depending on the result of the test
#' @examples
#' # No example yet
#' @export
check_adapted <- function(samples, unq_vals = 20) {
  # Only need to check uniqueness for one parameter
  first_par <- samples[1, , ]
  all(
    lapply(
      apply(first_par, 1, unique),
      length
    ) > unq_vals
  )
}
