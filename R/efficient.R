#' Check for efficient proposals if necessary
#'
#' Takes a mix proportion vector (3 x float) and the efficient proposal mu and
#' sigma. If efficient proposals are to be used (mix_proportion[3] > 0) then
#' test the efficient proposal values to see whether they are not null and
#' appropriate.
#'
#' @param efficient_mu The mu value for the efficient proposals
#' @param efficient_sig2 The sigma value for the efficient proposals
#' @param mix_proportion A vector of floats between 0 and 1 and summing to 1
#'   which give the proportion of particles to generate from the population
#'   level parameters, the individual random effects and the conditional
#'   parameters respectively
#'
#' @return nothing, stops operation on incorrect combination of parameters.
#' @keywords internal
check_efficient <- function(mix_proportion, efficient_mu, efficient_sig2) {
  if (mix_proportion[3] != 0) {
    if (is.null(efficient_mu) || is.null(efficient_sig2)) {
      stop(
        paste0(
          "Mu and sigma from efficient conditional ",
          "proposals must be provided for mix_proportion[3] > 0"
        )
      )
    }
  }
}


#' Create distribution parameters for efficient proposals
#'
#' From the existing samples, create a proposal distribution for drawing
#' efficient samples from.
#'
#' @param x The current pmwgs object
#'
#' @return A list containing the mu and sigma for the proposal distribution.
#' @keywords internal
create_efficient <- function(x) {
  proposal_means <- array(dim = c(x$n_pars, x$n_subjects))
  proposal_sigmas <- array(dim = c(x$n_pars, x$n_pars, x$n_subjects))
  for (s in 1:x$n_subjects) {
    cparms <- conditional_parms(
      s,
      extract_samples(x)
    )
    proposal_means[, s] <- cparms$cmeans
    proposal_sigmas[, , s] <- cparms$cvars
  }
  list(
    efficient_mu = proposal_means,
    efficient_sig2 = proposal_sigmas
  )
}


#' Obtain the efficent mu and sigma from the adaptation phase draws
#'
#' @param s current subject number
#' @param samples A list containing previous samples
#'
#' @return A list containing the conditional mean and variances for this subject
#' @keywords internal
conditional_parms <- function(s, samples) {
  tmudim <- dim(samples$theta_mu)
  n_par <- tmudim[1]
  n_iter <- tmudim[2]
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
  if (!isSymmetric(sigma_tilde))
    stop("covariance matrix for subject #", s, " not symmetric")
  if (any(eigen(sigma_tilde, TRUE, only.values = TRUE)$values < 1e-08))
    stop("covariance matrix for subject #", s, " not positive definite")
  condmvn <- condMVNorm::condMVN(
    mean = mu_tilde,
    sigma = sigma_tilde,
    dependent.ind = 1:n_par,
    given.ind = (n_par + 1):length(mu_tilde),
    X.given = c(
      samples$theta_mu[, n_iter],
      unwind(samples$theta_sig[, , n_iter])
    ),
    check.sigma = FALSE
  )
  list(cmeans = condmvn$condMean, cvars = condmvn$condVar)
}


#' Test that the sampler has successfully adapted
#'
#' @param pmwgs The full pmwgs object with all samples
#' @param n_unique The number of unique samples to look for in random effects
#'   for each subject.
#' @param i The number for the current iteration of the sampler
#'
#' @return A list containing a string representing successful/unsuccessful
#'   adaptation and an optional message. The string representing the success
#'   or failure can be one of c("success", "continue", "increase")
#'
#' @keywords internal
test_sampler_adapted <- function(pmwgs, n_unique, i) {
  fail_msg <- "values used in failed attempt to create proposal distribution"
  succ_msg <- "iterations before successful adaptation"
  if (i < n_unique) {
    return("continue")
  }
  test_samples <- extract_samples(pmwgs, stage = "adapt")
  if (check_adapted(test_samples$alpha, unq_vals = n_unique)) {
    attempt <- try({
      lapply(
        X = 1:pmwgs$n_subjects,
        FUN = conditional_parms,
        samples = test_samples
      )
    },
    silent = TRUE)
    if (inherits(attempt, "try-error")) {
      return(list("increase", paste("WARNING:", n_unique, fail_msg, "\n")))
    } else {
      return(list("success", paste("MESSAGE:", i, succ_msg, "\n")))
    }
  }
  return(list("continue"))
}


#' Check whether the adaptation phase has successfully completed
#'
#' @param samples The subject mean samples with which we are working
#' @param unq_vals The number of unique values for each subject
#'
#' @return A boolean TRUE or FALSE depending on the result of the test
#' @keywords internal
check_adapted <- function(samples, unq_vals = 20) {
  # Only need to check uniqueness for one parameter
  first_par <- samples[1, , ]
  # Split the matrix into a list of vectors by subject
  # Needed for the case where every sample is unique for all subjects
  first_par_list <- split(first_par, seq(NROW(first_par)))
  # Get unique pars (new accepted particles) and check length for
  # all subjects is greater than unq_vals
  all(lapply(lapply(first_par_list, unique), length) > unq_vals)
}
