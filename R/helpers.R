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
#' @param init list containing previous samples and more
#' @param ptm a single iteration version of the param_theta_mu variable
#' @param pts2 a single iter version of the param_theta_sigma2 variable
#' @param s current subject number
#' @param start_idx index specifying the starting sample to draw conditional distribution from
#' @param end_idx index specifying the end sample to draw coditional distribution from
#'
#' @return A list containing the conditional mean and variances for this subject
#' @examples
#' # No example yet
#' @export
conditional_parms <- function(init, ptm, pts2, s, start_idx, end_idx) {
  pts2_unwound <- apply(
    init$param_theta_sigma2[, , start_idx:end_idx],
    3,
    unwind
  )
  all_samples <- rbind(
    init$latent_theta_mu[, s, start_idx:end_idx],
    init$param_theta_mu[, start_idx:end_idx],
    pts2_unwound
  )
  mu_tilde <- apply(all_samples, 1, mean)
  sigma_tilde <- stats::var(t(all_samples))
  condmvn <- condMVNorm::condMVN(
    mean = mu_tilde,
    sigma = sigma_tilde,
    dependent.ind = 1:length(ptm),
    given.ind = (length(ptm) + 1):length(mu_tilde),
    X.given = c(ptm, unwind(pts2))
  )
  list(cmeans = condmvn$condMean, cvars = condmvn$condVar)
}


#' Create a new list for storage samples in the pmwgs object
#'
#' @param par_names The names of each parameter as a character vector
#' @param n_subjects The number of subjects for the subject mean storage.
#' @param iters The number of iterations to be pre-allocated
#'
#' @return A list containing the conditional mean and variances for this subject
#' @examples
#' # No example yet
#' @export
sample_store <- function(par_names, n_subjects, iters = 1) {
  n_pars <- length(par_names)
  list(
    subject_mean = array(
      NA,
      dim = c(n_pars, n_subjects, iters),
      dimnames = list(par_names, NULL, NULL)
    ),
    group_mean = array(
      NA,
      dim = c(n_pars, iters),
      dimnames = list(par_names, NULL)
    ),
    group_var = array(
      NA,
      dim = c(n_pars, n_pars, iters),
      dimnames = list(par_names, par_names, NULL)
    )
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
  if (anyNA(store$group_mean)) {
    last_ind <- which(
      apply(store$group_mean, 2, is.na),
      arr.ind = TRUE
    )[1, "col"]
  }
  else {
    last_ind <- ncol(store$group_mean)
  }

  list(
    gm = store$group_mean[, last_ind],
    gv = store$group_var[, , last_ind],
    sm = store$subject_mean[, , last_ind],
    gvi = store$last_group_var_inverse
  )
}
