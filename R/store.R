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
  stage_filter <- samples$stage %in% stage
  sampled_filter <- seq_along(samples$stage) <= samples$idx

  list(
    theta_mu = samples$theta_mu[, stage_filter & sampled_filter],
    theta_sig = samples$theta_sig[, , stage_filter & sampled_filter],
    alpha = samples$alpha[, , stage_filter & sampled_filter]
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
#' @param subject_ids The unique ID of each subjects as a character vector
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
    epsilon = array(
      NA_real_,
      dim = c(n_subjects, iters),
      dimnames = list(subject_ids, NULL)
    ),
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


#' Extend the main data store with empty space for new samples
#'
#' @param sampler The pmwgs object that we are adding the new samples to
#' @param n_samples The number of new samples to increase by
#' @param stage The name of the stage from which the new samples will be drawn
#'
#' @return The pmwgs object with the space for new samples added
#' @keywords internal
extend_sampler <- function(sampler, n_samples, stage) {
  old <- sampler$samples
  par_names <- sampler$par_names
  subject_ids <- sampler$subjects
  start <- old$idx + 1
  end <- old$idx + n_samples

  new_tmu <- array(NA_real_,
                   dim = dim(old$theta_mu) + c(0, n_samples),
                   dimnames = list(par_names, NULL))
  new_tmu[, - (start:end)] <- old$theta_mu
  sampler$samples$theta_mu <- new_tmu

  new_tsig <- array(NA_real_,
                    dim = dim(old$theta_sig) + c(0, 0, n_samples),
                    dimnames = list(par_names, par_names, NULL))
  new_tsig[, , - (start:end)] <- old$theta_sig
  sampler$samples$theta_sig <- new_tsig

  new_alph <- array(NA_real_,
                    dim = dim(old$alpha) + c(0, 0, n_samples),
                    dimnames = list(par_names, subject_ids, NULL))
  new_alph[, , - (start:end)] <- old$alpha
  sampler$samples$alpha <- new_alph

  new_epsilon <- array(NA_real_,
                   dim = dim(old$epsilon) + c(0, n_samples),
                   dimnames = list(subject_ids, NULL))
  new_epsilon[, - (start:end)] <- old$epsilon
  sampler$samples$epsilon <- new_epsilon

  new_sll <- array(NA_real_,
                   dim = dim(old$subj_ll) + c(0, n_samples),
                   dimnames = list(subject_ids, NULL))
  new_sll[, - (start:end)] <- old$subj_ll
  sampler$samples$subj_ll <- new_sll

  new_ahalf <- array(NA_real_,
                     dim = dim(old$a_half) + c(0, n_samples),
                     dimnames = list(par_names, NULL))
  new_ahalf[, - (start:end)] <- old$a_half
  sampler$samples$a_half <- new_ahalf

  sampler$samples$stage <- c(old$stage, rep(stage, n_samples))
  sampler
}


#' Trim the unneeded NA values from the end of the sampler
#'
#' @param sampler The pmwgs object that we are adding the new samples to
#'
#' @return The pmwgs object without NA values added during extend_sampler
#' @keywords internal
trim_na <- function(sampler) {
  idx <- sampler$samples$idx
  sampler$samples$theta_mu <- sampler$samples$theta_mu[, 1:idx]
  sampler$samples$theta_sig <- sampler$samples$theta_sig[, , 1:idx]
  sampler$samples$alpha <- sampler$samples$alpha[, , 1:idx]
  sampler$samples$subj_ll <- sampler$samples$subj_ll[, 1:idx]
  sampler$samples$a_half <- sampler$samples$a_half[, 1:idx]
  sampler$samples$stage <- sampler$samples$stage[1:idx]
  sampler$samples$epsilon <- sampler$samples$epsilon[, 1:idx]
  sampler
}


#' Relabel requested burn-in samples as adaptation
#'
#' Given a sampler object and a vector of sample indices, relabel the given
#' samples to be adaptation samples. This will allow them to be used in the
#' calculation of the conditional distribution for efficient sampling.
#'
#' @section Further information:
#'
#' This should not usually be needed, however if you have a model that is slow
#' to fit, and upon visual inspection and/or trace analysis you determine that
#' during burn-in the samples had already approached the posterior distribution
#' then you can use this function to re-label samples from that point onwards
#' to be classed as adaptation samples.
#'
#' This will allow them to be used in tests that check for the number of unique
#' samples, and in the building of the conditional distribution (which is used
#' for efficient sampling)
#'
#' If all old samples do not match `from` then an error will be raised.
#'
#' @param sampler The pmwgs object that we are relabelling samples from
#' @param indices The sample iterations from burn-in to relabel
#' @param from The stage you want to re-label from. Default is "burn"
#' @param to The stage you want to relabel to. Default is "adapt"
#'
#' @return The pmwgs object with re-labelled samples
#' @examples
#' new_pmwgs <- relabel_samples(sampled_forstmann, 17:21)
#' @export
relabel_samples <- function(sampler, indices, from = "burn", to = "adapt") {
  old_stage <- sampler$samples$stage
  if (!all(old_stage[indices] %in% from)) {
    stop(paste("Not all samples were from the", from, "stage"))
  }
  sampler$samples$stage[indices] <- to
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
#'         sampling process, i.e. one or more of "init", "burn", "adapt" or
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
#' par_estimates <- as_mcmc(sampled_forstmann)
#' par_estimates_sample_stage <- as_mcmc(sampled_forstmann, filter = "sample")
#' rand_eff <- as_mcmc(
#'   sampled_forstmann,
#'   selection = "alpha",
#'   filter = "sample"
#' )
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

#' Augment existing sampler object to have subject specific epsilon storage
#'
#' Older sampler object will be missing subject specific scaling parameter
#' (epsilon) storage, and running a stage with an updated pmwg will fail. To
#' fix this you can run the augment_sampler_epsilon function to fill the
#' appropriate array internals with NA values
#'
#' @param sampler The sampler object to augment
#'
#' @return A pmwgs sampler with epsilon array set internally
#'
#' @export
augment_sampler_epsilon <- function(sampler) {
  new_epsilon <- array(NA_real_,
                   dim = c(sampler$n_subjects, sampler$samples$idx),
                   dimnames = list(sampler$subjects, NULL))
  sampler$samples$epsilon <- new_epsilon
  sampler
}
