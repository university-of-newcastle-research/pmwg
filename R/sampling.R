#' Run a stage of the PMwG sampler
#'
#' Run one of burnin, adaptation or sampling phase from the PMwG
#' sampler. Each stage involves slightly different processes, so for the
#' full PMwG sampling we need to run this three times.
#'
#' The \strong{burnin} stage by default selects 50% of particles from the model
#' parameter sample (selected through a Gibbs step) and 50% of particles from
#' the previous random effect of each participant. It assesses each particle
#' with the log-likelihood function and samples from all particles weighted by
#' their log-likelihood.
#'
#' The \strong{adaptation} stage selects and assesses particle in the same was
#' as burnin, however on each iteration it also checks whether each participant
#' has enough unique random effect samples to attempt to create a conditional
#' distribution for efficient sampling in the next stage. If the attempt at
#' creating a conditional distribution fails, then the number of unique samples
#' is increased and sampling continues. If the attempt succeeds then the stage
#' is finished early.
#'
#' The \strong{final} stage (sampling) by default samples predominantly from the
#' conditional distribution created at the end of adaptation. This is more
#' efficient and allows the number of particles to be reduced whilst still
#' getting a high enough acceptance rate of new samples.
#'
#' Once complete each stage will return a sampler object with the new samples
#' stored within it.
#'
#' @param pmwgs A Particle Metropolis within Gibbs sampler which has been set
#'   up and initialised
#' @param stage The sampling stage to run. Must be one of \code{'burn'},
#'   \code{'adapt'} or \code{'sample'}.
#' @param iter The number of iterations to run for the sampler. For
#'   \code{'burn'} and \code{'sample'} all iterations will run. However for
#'   \code{'adapt'} if all subjects have enough unique samples to create the
#'   conditional distribution then the stage will finish early.
#' @param particles The default here is 1000 particles to be generated for each
#'   iteration, however during the sample phase this should be reduced.
#' @param display_progress Display a progress bar during sampling.
#' @param n_cores Set to more than 1 to use \code{mclapply}. Setting
#'   \code{n_cores} greater than 1 is only permitted on Linux and Mac OS X
#'   machines.
#' @param ... Further arguments used to fine tune the sampling process. Accepted
#'   additional arguments include:
#'   \itemize{
#'     \item \code{epsilon} which should be a value between 0 and 1.
#'       \strong{epsilon} controls the extent to which the covariance matrix is
#'       scaled when generating particles from the previous random effect.
#'     \item \code{mix} controls the mixture of different sources for particles.
#'       The function \code{\link{numbers_from_proportion}} is passed this value
#'       and includes extra details on what is accepted.
#'     \item \code{n_unique} is a number representing the number of unique
#'       samples to check for on each iteration of the sampler. Only used during
#'       the \code{'adapt'} stage.
#'   }
#'
#' @return A pmwgs object with the newly generated samples in place.
#' @examples
#' # No example yet
#' @export
run_stage <- function(pmwgs,
                      stage,
                      iter = 1000,
                      particles = 1000,
                      display_progress = TRUE,
                      n_cores = 1,
                      ...) {
  # Check the passed arguments and return cleaned - extract relevant arguments
  .args <- as.list(match.call()[-1])
  clean_args <- do.call(check_run_stage_args, .args)
  stage <- clean_args$stage
  if (stage == "adapt") {
    .n_unique <- clean_args$.n_unique
    n_unique <- clean_args$n_unique
  }
  apply_fn <- clean_args$apply_fn

  # Display stage to screen
  msgs <- list(
    burn = "Phase 1: Burn in\n", adapt = "Phase 2: Adaptation\n",
    sample = "Phase 3: Sampling\n"
  )
  cat(msgs[[stage]])

  # Build new sample storage
  stage_samples <- sample_store(
    pmwgs$par_names, pmwgs$subjects,
    iters = iter, stage = stage
  )
  # create progress bar
  if (display_progress) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }

  # Main iteration loop
  for (i in 1:iter) {
    if (display_progress) {
      update_progress_bar(pb, i, extra = mean(accept_rate(stage_samples)))
    }

    if (i == 1) store <- pmwgs$samples else store <- stage_samples
    tryCatch(
      pars <- gibbs_step(store, pmwgs),
      error = function(err_cond) {
        gibbs_step_err(pmwgs, store)
      }
    )

    # Sample new particles for random effects.
    # Send new_sample the "index" of the subject id - not subject id itself.
    pmwgs_args <- list(
      X = 1:pmwgs$n_subjects,
      FUN = new_sample,
      data = pmwgs$data,
      num_particles = particles,
      parameters = pars,
      likelihood_func = pmwgs$ll_func,
      subjects = pmwgs$subjects
    )
    fn_args <- c(pmwgs_args, clean_args$sample_args)
    tmp <- do.call(apply_fn, fn_args)

    ll <- unlist(lapply(tmp, attr, "ll"))
    alpha <- array(unlist(tmp), dim = dim(pars$alpha))

    # Store results locally.
    stage_samples$theta_mu[, i] <- pars$tmu
    stage_samples$theta_sig[, , i] <- pars$tsig
    stage_samples$last_theta_sig_inv <- pars$tsinv
    stage_samples$alpha[, , i] <- alpha
    stage_samples$idx <- i
    stage_samples$subj_ll[, i] <- ll
    stage_samples$a_half[, i] <- pars$a_half

    if (stage == "adapt") {
      res <- test_sampler_adapted(stage_samples, pmwgs, n_unique, i)
      if (res == "success") {
        break
      } else if (res == "increase") {
        n_unique <- n_unique + .n_unique
      }
    }

  }
  if (display_progress) close(pb)
  if (stage == "adapt") {
    if (i == iter) {
      message(paste(
        "Particle Metropolis within Gibbs Sampler did not",
        "finish adaptation phase early (all", i, "iterations were",
        "run).\nYou should examine your samples and perhaps start",
        "a longer adaptation run."
      ))
    }
  }
  update_sampler(pmwgs, stage_samples)
}


#' Generate particles and select one to be the new sample
#'
#' Generate a new sample for a particular subject given their data and the
#' new model parameter estimates. This should not be called directly, rather it
#' is used internally to run_stage.
#'
#' The function that controls the generation of new samples for the Particle
#' Metropolis within Gibbs sampler. It generates samples from either the initial
#' proposal or from the last iteration of the sampler. This function should not
#' usually need to be called, as the \code{run_stage} function uses this
#' internally.
#'
#' The way it selects a new sample is by generating proposal particles from up
#' to three different distributions (according to a mixing proportion).
#'
#' The first distribution is based on the current model parameter sample values.
#' The second distribution is based on the last random effects for the subject.
#' The third distribution is only used in the final sampling phase and is based
#' on the conditional distribution built from accepted particles in the adapt
#' phase of the sampler.
#'
#' @param s A number - the index of the subject. For \code{s == 1} The first
#'   subject ID from the \code{data} subject column will be selected. For
#'   \code{s == 2} the second unique value for subject id will be used.
#' @param data A data.frame (or similar object) which contains the data against
#'   which the particles are assessed. The only strict requirement is that
#'   it contains a subject column named as such to allow for the splitting
#'   of the data by unique subject id. The provided log likelihood function
#'   is the only other contact with the data.
#' @param parameters A list containing:
#'   \describe{
#'     \item{\code{tmu}}{The vector of means for the multivariate normal}
#'     \item{\code{tsig}}{A covariate matrix for the multivariate normal}
#'     \item{\code{alpha}}{An array of individual subject random effects}
#'   }
#' @inheritParams numbers_from_proportion
#' @inheritParams check_efficient
#' @param likelihood_func A likelihood function for calculating log likelihood
#'   of samples. Usually provided internally in \code{run_stage} from the pmwgs
#'   object.
#' @param epsilon A scaling factor to reduce the variance on the distribution
#'   based on subject random effects when generating particles.
#' @param subjects A list of unique subject ids in the order they appear in
#'   the data.frame
#'
#' @return A single sample from the new proposals
#' @keywords internal
new_sample <- function(s, data, num_particles, parameters,
                       efficient_mu = NULL, efficient_sig2 = NULL,
                       mix_proportion = c(0.5, 0.5, 0.0),
                       likelihood_func = NULL,
                       epsilon = 1, subjects = NULL) {
  # Check for efficient proposal values if necessary
  check_efficient(mix_proportion, efficient_mu, efficient_sig2)
  e_mu <- efficient_mu[, s]
  e_sig2 <- efficient_sig2[, , s]
  mu <- parameters$tmu
  sig2 <- parameters$tsig
  subj_mu <- parameters$alpha[, s]
  if (is.null(likelihood_func)) stop("likelihood_func is a required argument")

  # Create proposals for new particles
  proposals <- gen_particles(
    num_particles, mu, sig2, subj_mu,
    mix_proportion = mix_proportion,
    prop_mu = e_mu,
    prop_sig2 = e_sig2,
    epsilon = epsilon
  )
  # Put the current particle in slot 1.
  proposals[1, ] <- subj_mu

  # Density of data given random effects proposal.
  lw <- apply(
    proposals,
    1,
    likelihood_func,
    data = data[data$subject == subjects[s], ]
  )
  # Density of random effects proposal given population-level distribution.
  lp <- mvtnorm::dmvnorm(x = proposals, mean = mu, sigma = sig2, log = TRUE)
  # Density of proposals given proposal distribution.
  prop_density <- mvtnorm::dmvnorm(
    x = proposals,
    mean = subj_mu,
    sigma = sig2 * (epsilon^2)
  )
  # Density of efficient proposals
  if (mix_proportion[3] != 0) {
    eff_density <- mvtnorm::dmvnorm(
      x = proposals,
      mean = e_mu,
      sigma = e_sig2
    )
  } else {
    eff_density <- 0
  }

  lm <- log(mix_proportion[1] * exp(lp) +
    (mix_proportion[2] * prop_density) +
    (mix_proportion[3] * eff_density))
  # log of importance weights.
  l <- lw + lp - lm
  weights <- exp(l - max(l))
  idx <- sample(x = num_particles, size = 1, prob = weights)
  winner <- proposals[idx, ]
  attr(winner, "ll") <- lw[idx]
  winner
}


#' Generate proposal particles
#'
#' Generates particles for the \code{new_sample} function
#'
#' Generate particles for a particular subject from a mix of population level
#' (hierarchical) distribution, from the particles (containing individual level
#' distribution) and/or from the conditional on accepted individual level
#' particles, a more efficient proposal method.
#'
#' This function is used in burnin, adaptation and sampling using various
#' combinations of the arguments.
#'
#' @param mu A vector of means for the multivariate normal
#' @param sig2 A covariate matrix for the multivariate normal
#' @param particle A particle (re proposals for latent variables)
#' @inheritParams numbers_from_proportion
#' @param epsilon Reduce the variance for the individual level samples by this
#'   factor
#'
#' @return The new proposals
#' @keywords internal
gen_particles <- function(num_particles,
                          mu,
                          sig2,
                          particle,
                          ...,
                          mix_proportion = c(0.5, 0.5, 0.0),
                          prop_mu = NULL,
                          prop_sig2 = NULL,
                          epsilon = 1) {
  particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
  # Generate proposal particles
  pop_particles <- particle_draws(particle_numbers[1], mu, sig2)
  ind_particles <- particle_draws(
    particle_numbers[2],
    particle,
    sig2 * epsilon^2
  )
  eff_particles <- particle_draws(particle_numbers[3], prop_mu, prop_sig2)
  particles <- rbind(pop_particles, ind_particles, eff_particles)
  colnames(particles) <- names(mu) # stripped otherwise.
  particles
}


#' Check and normalise the number of each particle type from the mix_proportion
#'
#' Takes a mix proportion vector (3 x float) and a number of particles to
#' generate and returns a vector containing the number of each particle type to
#' generate.
#'
#' @param mix_proportion A vector of floats between 0 and 1 and summing to 1
#'   which give the proportion of particles to generate from the population
#'   level parameters, the individual random effects and the conditional
#'   parameters respectively
#' @param num_particles The total number of particles to generate using a
#'   combination of the three methods.
#'
#' @return The wound vector as a matrix
#' @keywords internal
numbers_from_proportion <- function(mix_proportion, num_particles = 1000) {
  numbers <- stats::rbinom(n = 2, size = num_particles, prob = mix_proportion)
  if (mix_proportion[3] == 0) {
    numbers[3] <- 0
    numbers[2] <- num_particles - numbers[1]
  } else {
    numbers[3] <- num_particles - sum(numbers)
  }
  numbers
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
#' @keywords internal
particle_draws <- function(n, mu, covar) {
  if (n <= 0) {
    return(NULL)
  }
  mvtnorm::rmvnorm(n, mean = mu, sigma = covar)
}


#' Test the arguments to the run_stage function for correctness
#'
#' Takes the arguments to run_stage and checks them for completeness and
#' correctness. Returns a list of cleaned/checked arguments to the caller.
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
  sargs <- list()
  # Check pmwgs object is correct type and has been initialised
  if (!is.pmwgs(pmwgs)) {
    stop("`run_stage` function Requires an object of type <pmwgs>")
  }
  if (!pmwgs$init) stop("pmwgs object has not been initialised")

  # Test stage argument
  valid_stages <- stages[2:length(stages)]
  stage_err <- paste(c("Argument `stage` should be one of -", valid_stages),
                     collapse = " ",
                     sep = ",")
  tryCatch(
    stage <- match.arg(stage, valid_stages),
    error = function(err_cond) {
      stop(stage_err)
    }
  )

  acceptable_extras <- c("n_unique", "epsilon", "mix")
  dots <- list(...)
  for (argname in names(dots)) {
    if (!argname %in% acceptable_extras) {
      stop(paste("Unexpected argument", argname, "passed to run_stage"))
    }
  }
  # Extract n_unique argument
  if (stage == "adapt") {
    adapt_args <- list()
    if (is.null(dots$n_unique)) {
      adapt_args$.n_unique <- 20
    } else {
      adapt_args$.n_unique <- dots$n_unique
      dots$n_unique <- NULL
    }
    adapt_args$n_unique <- adapt_args$.n_unique
  } else {
    adapt_args <- list()
    if (!is.null(dots$n_unique)) {
      dots$n_unique <- NULL
      warning("Argument `n_unique` unused for any stage other than adapt")
    }
  }

  # Set a default value for epsilon if it does not exist
  if (is.null(dots$epsilon)) {
    sargs$epsilon <- ifelse(pmwgs$n_pars > 15,
      0.1,
      ifelse(pmwgs$n_pars > 10, 0.3, 0.5)
    )
    message(paste("Epsilon has been set automatically to:", sargs$epsilon))
    message("This is based on the number of parameters")
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
    sargs <- c(sargs, prop_args)
  }

  # Set default values for the mix_proportion parameter if not passed in as arg,
  # and perform checks on its values/length
  if (is.null(dots$mix)) {
    if (stage == "sample") {
      sargs$mix <- c(0.1, 0.2, 0.7)
    } else {
      sargs$mix <- c(0.5, 0.5, 0.0)
    }
  } else {
    sargs$mix <- dots$mix
  }
  if (!isTRUE(all.equal(sum(sargs$mix), 1))) {
    stop("The elements of the `mix` proportion vector must sum to 1")
  }
  if (length(sargs$mix) != 3) {
    stop("`mix` proportion vector must have three elements which sum to 1")
  }

  apply_fn <- lapply
  if (n_cores > 1) {
    if (Sys.info()[["sysname"]] == "Windows") {
      stop("`n_cores` cannot be greater than 1 on Windows systems.")
    }
    apply_fn <- parallel::mclapply
    sargs$mc.cores <- n_cores #nolint
  }
  c(adapt_args, list(stage = stage, sample_args = sargs, apply_fn = apply_fn))
}


#' Return the acceptance rate for all subjects
#'
#' @param store The samples store (containing random effects) with which we are
#'   working
#'
#' @return A vector with the acceptance rate for each subject
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
