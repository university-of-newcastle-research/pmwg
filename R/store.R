#' Extract relevant samples from the list for conditional dist calc
#'
#' From the sampler, extract relevant samples for the creation of
#' the proposal distribution.
#'
#' @param sampler The pmwgs object containing the samples
#' @param stage The stage, or list of stages from which you want the samples
#'
#' @return A list containing only appropriate samples (non init/burnin samples)
#' @keywords internal
extract_samples <- function(sampler, stage = c("adapt", "sample")) {
  samples <- sampler$samples
  sample_filter <- samples$stage %in% stage
  list(
    theta_mu = samples$theta_mu[, sample_filter],
    theta_sig = samples$theta_sig[, , sample_filter],
    alpha = samples$alpha[, , sample_filter]
  )
}


#' Unwinds variance matrix to a vector
#'
#' Takes a variance matrix and unwind to a vector via Cholesky decomposition
#' then take the log of the diagonal.
#'
#' @param var_matrix A variance matrix
#'
#' @return The unwound matrix as a vector
#' @keywords internal
unwind <- function(var_matrix, ...) {
  y <- t(chol(var_matrix))
  diag(y) <- log(diag(y))
  y[lower.tri(y, diag = TRUE)]
}

#' Winds a variance vector back to a vector
#'
#' The reverse of the unwind function, takes a variance vector and windows back
#' into matrix
#'
#' @param var_vector A variance vector
#'
#' @return The wound vector as a matrix
#' @keywords internal
wind <- function(var_vector, ...) {
  n <- sqrt(2 * length(var_vector) + 0.25) - 0.5 ## Dim of matrix.
  if ((n * n + n) != (2 * length(var_vector))) stop("Wrong sizes in unwind.")
  out <- array(0, dim = c(n, n))
  out[lower.tri(out, diag = TRUE)] <- var_vector
  diag(out) <- exp(diag(out))
  out %*% t(out)
}


#' Create a new list for storage samples in the pmwgs object
#'
#' @param par_names The names of each parameter as a character vector
#' @param subjects_ids The unique ID of each subjects as a character vector
#' @param iters The number of iterations to be pre-allocated
#' @param stage The stage for which the samples will be created. Should be one
#'   of \code{c("init", "burn", "adapt", "sample")}
#'
#' @return A list containing the conditional mean and variances for this subject
#' @keywords internal
sample_store <- function(par_names, subject_ids, iters = 1, stage = "init") {
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  list(
    alpha = array(
      NA_real_,
      dim = c(n_pars, n_subjects, iters),
      dimnames = list(par_names, subject_ids, NULL)
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
    stage = array(stage, iters),
    subj_ll = array(
      NA_real_,
      dim = c(n_subjects, iters),
      dimnames = list(subject_ids, NULL)
    ),
    a_half = array(
      NA_real_,
      dim = c(n_pars, iters),
      dimnames = list(par_names, NULL)
    )
  )
}


#' Update the main data store with the results of the last stage
#'
#' @param sampler The pmwgs object that we are adding the new samples to
#' @param store The sample storage stage just run
#'
#' @return The pmwgs object with the new samples concatenated to the old
#' @keywords internal
update_sampler <- function(sampler, store) {
  old_tmu <- sampler$samples$theta_mu
  old_tsig <- sampler$samples$theta_sig
  old_alpha <- sampler$samples$alpha
  old_stage <- sampler$samples$stage
  old_sll <- sampler$samples$subj_ll
  old_a_half <- sampler$samples$a_half
  li <- store$idx
  par_names <- sampler$par_names
  subject_ids <- sampler$subjects

  sampler$samples$theta_mu <- array(
    c(old_tmu, store$theta_mu[, 1:li]),
    dim = dim(old_tmu) + c(0, li),
    dimnames = list(par_names, NULL)
  )
  sampler$samples$theta_sig <- array(
    c(old_tsig, store$theta_sig[, , 1:li]),
    dim = dim(old_tsig) + c(0, 0, li),
    dimnames = list(par_names, par_names, NULL)
  )
  sampler$samples$alpha <- array(
    c(old_alpha, store$alpha[, , 1:li]),
    dim = dim(old_alpha) + c(0, 0, li),
    dimnames = list(par_names, subject_ids, NULL)
  )
  sampler$samples$idx <- ncol(sampler$samples$theta_mu)
  sampler$samples$last_theta_sig_inv <- store$last_theta_sig_inv
  sampler$samples$stage <- c(old_stage, store$stage[1:li])
  sampler$samples$subj_ll <- array(
    c(old_sll, store$subj_ll[, 1:li]),
    dim = dim(old_sll) + c(0, li),
    dimnames = list(subject_ids, NULL)
  )
  sampler$samples$a_half <- array(
    c(old_a_half, store$a_half[, 1:li]),
    dim = dim(old_a_half) + c(0, li),
    dimnames = list(par_names, NULL)
  )
  sampler
}


#' Create a list with the last samples in the pmwgs object
#'
#' @param store The list containing samples from which to grab the last.
#'
#' @return A list containing the last sample of group mean and variance and
#'   subject means.
#' @keywords internal
last_sample <- function(store) {
  list(
    tmu = store$theta_mu[, store$idx],
    tsig = store$theta_sig[, , store$idx],
    alpha = store$alpha[, , store$idx],
    tsinv = store$last_theta_sig_inv,
    a_half = store$a_half[, store$idx]
  )
}


#' Return a CODA mcmc object with the required samples
#'
#' Given a sampler object and a specification of the samples required, return
#' either an individual coda mcmc object, or a list of mcmc objects.
#'
#' @section Selecting sample types:
#'
#' The values that can be chosen for the \code{selection} argument can be one
#' of the following list:
#' \describe{
#'   \item{\code{"theta_mu"}}{the model parameter estimate samples}
#'   \item{\code{"theta_sig"}}{the covariance matrix estimates, returns a list
#'     of mcmc objects, one for each model parameter.}
#'   \item{\code{"alpha"}}{the random effect estimates, returns a list of mcmc
#'     objects, one for each subject.}
#' }
#' The default value for \code{selection} is \code{"theta_mu"}
#'
#' @section Filtering samples:
#'
#' The \code{filter} argument can take one of two forms:
#' \itemize{
#'   \item An integer vector, usually a sequence of integers, that must fall
#'         within the range 1:end.
#'   \item A character vector, where each element corresponds to a stage of the
#'         sampling process, ie one or more of "init", "burn", "adapt" or
#'         "sample".
#' }
#' The default value for \code{filter} is all stages.
#'
#' @param sampler The pmwgs object containing samples to extract.
#' @param selection The selection of sample types to return.
#' @param filter A filter that defines which stage to draw samples from.
#'
#' @return An mcmc object or list containing the selected samples.
#' @examples
#' # No example yet
#' @export
as_mcmc <- function(sampler, selection = "theta_mu", filter = stages) {
  if (all(filter %in% stages)) {
    filter <- which(sampler$samples$stage %in% filter)
  } else if (!all(filter %in% 1:sampler$samples$idx)) {
    stop("filter is not a vector of stage names, or integer vector of indices")
  }

  if (selection == "theta_mu") {
    return(coda::mcmc(t(sampler$samples$theta_mu[, filter])))
  } else if (selection == "theta_sig") {
    tsig <- sampler$samples$theta_sig[, , filter]
    return(stats::setNames(lapply(
      seq(dim(tsig)[1]),
      function(x) {
        coda::mcmc(t(tsig[x, , ]))
      }
    ), sampler$par_names))
  } else if (selection == "alpha") {
    alpha <- sampler$samples$alpha[, , filter]
    return(stats::setNames(lapply(
      seq(dim(alpha)[2]),
      function(x) {
        coda::mcmc(t(alpha[, x, ]))
      }
    ), sampler$subjects))
  }
  stop("Argument `selection` should be one of theta_mu, theta_sig, alpha")
}
