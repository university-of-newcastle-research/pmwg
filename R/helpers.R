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
#' The reverse of the unwind function, takes a variance vector and windows back
#' into matrix
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


#' Test the arguments to the run_stage function for correctness
#'
#' Takes the arguments to run_stage and checks them for completeness and
#' correctness. Uses parent.frame to edit run_stage env directly
#'
#' @inheritParams run_stage
#'
#' @keywords internal
check_run_stage_args <- function(pmwgs,
                                 stage,
                                 iter = 1000,
                                 particles = 1000,
                                 display_progress = TRUE,
                                 n_cores = 1,
                                 ...) {
  # Two lists of arguments used in different places to return
  run_args <- list()
  sample_args <- list()
  # Check pmwgs object is correct type and has been initialised
  pmwgs <- pmwgs
  if (!is.pmwgs(pmwgs)) {
    stop("`run_stage` function Requires an object of type <pmwgs>")
  }
  if (!pmwgs$init) stop("pmwgs object has not been initialised")

  # Test stage argument
  tryCatch(
    run_args$stage <- match.arg(stage, c("burn", "adapt", "sample")),
    error = function(err_cond) {
      stop("Argument `stage` should be one of 'burn', 'adapt' or 'sample'")
    }
  )

  dots <- list(...)
  # Extract n_unique argument
  if (stage == "adapt") {
    if (is.null(dots$n_unique)) {
      run_args$.n_unique <- 20
    } else {
      run_args$.n_unique <- dots$n_unique
      dots$n_unique <- NULL
    }
    run_args$n_unique <- .n_unique
  } else {
    if (!is.null(dots$n_unique)) {
      dots$n_unique <- NULL
      warning("Argument `n_unique` unused for any stage other than adapt")
    }
  }

  # Set a default value for epsilon if it does not exist
  if (is.null(dots$epsilon)) {
    sample_args$epsilon <- ifelse(pmwgs$n_pars > 15,
      0.1,
      ifelse(pmwgs$n_pars > 10, 0.3, 0.5)
    )
  }

  # Create efficient proposal distribution if we are in sampling phase
  if (stage == "sample") {
    tryCatch(
      prop_args <- try(create_efficient(pmwgs)),
      error = function(err_cond) {
        outfile <- tempfile(
          pattern = "PMwG_err_",
          tmpdir = ".",
          fileext = ".RData"
        )
        msg <- paste(
          "An error was detected whilst creating conditional",
          "distribution.\n",
          "Saving current state of environment in file:",
          outfile
        )
        save.image(outfile)
        stop(msg)
      }
    )
    run_args <- c(run_args, prop_args)
  }

  # Set default values for the mix_ratio parameter if not passed in as arg, and
  # perform checks on its values/length
  if (is.null(dots$mix)) {
    if (stage == "sample") {
      sample_args$mix <- c(0.1, 0.2, 0.7)
    } else {
      sample_args$mix <- c(0.5, 0.5, 0.0)
    }
  } else {
    sample_args$mix <- dots$mix
  }
  if (!isTRUE(all.equal(sum(sample_args$mix), 1))) {
    stop("The elements of the `mix` ratio vector must sum to 1")
  }
  if (length(sample_args$mix) != 3) {
    stop("`mix` ratio vector must have three elements which sum to 1")
  }

  apply_fn <- lapply
  if (n_cores > 1) {
    if (Sys.info()[["sysname"]] == "Windows") {
      stop("`n_cores` cannot be greater than 1 on Windows systems.")
    }
    apply_fn <- parallel::mclapply
    sample_args$mc.cores <- n_cores #nolint
  }
  list(run_args = run_args, sample_args = sample_args, apply_fn = apply_fn)
}


#' Check and normalise the number of each particle type from the mix_ratio
#'
#' Takes a mix ratio vector (3 x float) and a number of particles to generate
#' and returns a vector containing the number of each particle type to generate
#'
#' @param mix_ratio A vector of floats betwen 0 and 1 and summing to 1 which
#'   give the ratio of particles to generate from the population level
#'   parameters, the individual random effects and the conditional parameters
#'   repectively
#' @param num_particles The total number of particles to generate using a
#'   combination of the three methods.
#'
#' @return The wound vector as a matrix
#' @examples
#' psamplers:::numbers_from_ratio(c(0.1, 0.3, 0.6))
#' @keywords internal
numbers_from_ratio <- function(mix_ratio, num_particles = 1000) {
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
#' @param mix_ratio A vector of floats betwen 0 and 1 and summing to 1 which
#'   give the ratio of particles to generate from the population level
#'   parameters, the individual random effects and the conditional parameters
#'   repectively
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
#' From the sampler, extract relevant samples for the creation of
#' the proposal distribution.
#'
#' @param sampler The pmwgs object containing the samples
#' @param stage The stage, or list of stages from which you want the samples
#'
#' @return A list containing only appopriate samples (non init/burnin samples)
#' @examples
#' # No example yet
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


#' Generate a cloud of particles from a multivariate normal distribution
#'
#' Takes the mean and variance for a multivariate normal distribution, as well
#' as the number of particles to generate and return random draws from the
#' multivariate normal if the numbers of particles is > 0, otherwise return
#' NULL. At least one of mean or sigma must be provided.
#'
#' @param n number of observations
#' @param mu mean vector
#' @param covar covariance matrix
#'
#' @return If n > 0 returns n draws from the multivariate normal with mean and
#'   sigma, otherwise returns NULL
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
#' @param s current subject number
#' @param samples A list containing previous samples
#'
#' @return A list containing the conditional mean and variances for this subject
#' @examples
#' # No example yet
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
  condmvn <- condMVNorm::condMVN(
    mean = mu_tilde,
    sigma = sigma_tilde,
    dependent.ind = 1:n_par,
    given.ind = (n_par + 1):length(mu_tilde),
    # GC: Note, not sure what is happening here:v (Was ptm/pts2 now last sample)
    X.given = c(
      samples$theta_mu[, n_iter],
      unwind(samples$theta_sig[, , n_iter])
    )
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
#' @keywords internal
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
    stage = array(stage, iters),
    subj_ll = array(
      NA_real_,
      dim = c(n_subjects, iters),
      dimnames = list(NULL, NULL)
    ),
    a_half = array(
      NA_real_,
      dim = c(n_pars, iters),
      dimnames = list(par_names, NULL)
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


#' Update the main data store with the results of the last stage
#'
#' @param sampler The pmwgs object that we are adding the new samples to
#' @param store The sample storage stage just run
#'
#' @return The pmwgs object with the new samples concatenated to the old
#' @examples
#' # No example yet
#' @keywords internal
update_sampler <- function(sampler, store) {
  old_tmu <- sampler$samples$theta_mu
  old_tsig <- sampler$samples$theta_sig
  old_alpha <- sampler$samples$alpha
  old_stage <- sampler$samples$stage
  old_sll <- sampler$samples$subj_ll
  old_a_half <- sampler$samples$a_half
  li <- store$idx

  sampler$samples$theta_mu <- array(c(old_tmu, store$theta_mu[, 1:li]),
                                      dim = dim(old_tmu) + c(0, li))
  sampler$samples$theta_sig <- array(c(old_tsig, store$theta_sig[, , 1:li]),
                                     dim = dim(old_tsig) + c(0, 0, li))
  sampler$samples$alpha <- array(c(old_alpha, store$alpha[, , 1:li]),
                                        dim = dim(old_alpha) + c(0, 0, li))
  sampler$samples$idx <- ncol(sampler$samples$theta_mu)
  sampler$samples$last_theta_sig_inv <- store$last_theta_sig_inv
  sampler$samples$stage <- c(old_stage, store$stage[1:li])
  sampler$samples$subj_ll <- array(c(old_sll, store$subj_ll[, 1:li]),
                                   dim = dim(old_sll) + c(0, li))
  sampler$samples$a_half <- array(c(old_a_half, store$a_half[, 1:li]),
                                      dim = dim(old_a_half) + c(0, li))
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
#' @keywords internal
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


#' Return the acceptance rate for all subjects
#'
#' @param store The samples store (containing random effects) with which we are
#'   working
#'
#' @return A vector with the acceptance rate for each subject
#' @examples
#' # No example yet
#' @keywords internal
accept_rate <- function(store) {
  if (is.null(store$idx) || store$idx < 3) {
    return(array(0, dim(store$alpha)[2]))
  }
  vals <- store$alpha[1, , 1:store$idx]
  apply(
    apply(vals, 1, diff) != 0, # If diff != 0
    2,
    mean
  )
}


#' An altered version of the utils:txtProgressBar that shows acceptance rate
#'
#' @param min The minimum of the value being updated for the progress bar
#' @param max The maximum of the value being updated for the progress bar
#'
#' @return A structure matching the structure of a txtProgresBar with additional
#'   info
#' @keywords internal
accept_progress_bar <- function(min = 0, max = 1) {
  .val <- 0
  .killed <- FALSE
  .nb <- 0L
  .pc <- -1L # This ensures the initial value is displayed
  .ex <- 0
  nw <- nchar("=", "w")
  width <- trunc(getOption("width") - 22L / nw)
  if (max <= min) stop("must have 'max' > 'min'")

  up <- function(value, extra = 0) {
    if (!is.finite(value) || value < min || value > max) {
      return()
    }
    .val <<- value
    nb <- round(width * (value - min) / (max - min))
    pc <- round(100 * (value - min) / (max - min))
    extra <- round(100 * extra)
    if (nb == .nb && pc == .pc && .ex == extra) {
      return()
    }
    cat(paste0("\r  |", strrep(" ", nw * width + 6)))
    cat(paste(c(
      "\r  |",
      rep.int("=", nb),
      rep.int(" ", nw * (width - nb)),
      sprintf("| %3d%%", pc),
      sprintf(" | Acc(%3d%%)", extra)
    ), collapse = ""))
    utils::flush.console()
    .nb <<- nb
    .pc <<- pc
    .ex <<- extra
  }

  get_value <- function() .val
  kill <- function() {
    if (!.killed) {
      cat("\n")
      utils::flush.console()
      .killed <<- TRUE
    }
  }
  up(0) # will check if in range

  structure(list(getVal = get_value, up = up, kill = kill),
            class = c("accept_progress_bar", "txtProgressBar"))
}

#' A function that updates the accept_progress_bar with progress and accept rate
#'
#' @param pb The progress bar object
#' @param value The value to set the bar width to
#' @param extra A value that represents the number of accepted particles
#'
#' @return The old value that was present before updating
#'
#' @keywords internal
update_progress_bar <- function(pb, value, extra = 0) {
  if (!inherits(pb, "txtProgressBar")) {
    stop(gettextf(
      "'pb' is not from class %s",
      dQuote("txtProgressBar")
    ),
    domain = NA
    )
  }
  oldval <- pb$getVal()
  pb$up(value, extra)
  invisible(oldval)
}

#' Error handler for the new_group_pars call
#'
#' @param err_cond The samples store (containing random effects) with which we are
#'   working
#'
#' @return A vector with the acceptance rate for each subject
#' @examples
#' # No example yet
#' @keywords internal
new_group_pars_err <- function(pmwgs, store) {
  store_tmp <- tempfile(
    pattern = "pmwg_stage_samples_",
    tmpdir = ".",
    fileext = ".RDS"
  )
  sampler_tmp <- tempfile(
    pattern = "pmwg_obj_",
    tmpdir = ".",
    fileext = ".RDS"
  )
  message("Error while generating new group level parameters")
  message("Saving current state of pmwgs object: ", sampler_tmp)
  saveRDS(pmwgs, file = sampler_tmp)
  message("Saving current state of stage sample storage", store_tmp)
  saveRDS(store, file = store_tmp)
  stop("Stopping execution")
}
