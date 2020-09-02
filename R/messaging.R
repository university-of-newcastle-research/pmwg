
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


#' Error handler for the gibbs_step call
#'
#' If an error was detected when generating new values in Gibbs step this
#' function is called to generate the error message and save the state of the
#' samples at that moment to help with debugging.
#'
#' @param pmwgs The pmwgs object for the current run.
#' @param store The samples store (containing random effects) with which we
#'   are working in the current stage.
#'
#' @keywords internal
gibbs_step_err <- function(pmwgs, store) {
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
