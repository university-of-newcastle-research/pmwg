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
  if ( (n * n + n) != (2 * length(var_vector)) ) stop("Wrong sizes in unwind.")
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
