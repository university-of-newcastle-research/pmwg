#' An altered version of the utils:txtProgressBar that shows acceptance rate
#'
#' The progress bar displays several elements, the progress visually as a bar
#' being filled and the percentage complete as per the standard
#' utils::txtProgressBar and additionally the average across subjects of the
#' rate of accepting newly generated particles.
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
  component <- list(
    pchar = "=",
    prog_start = " |",
    prog_end = "| ",
    percent = "%3d%%",
    acc_sep = " | ",
    acc_msg = "New(%3d%%)"
  )
  width <- nchar(lapply(component, gettextf, 100), "w")
  width <- split(unname(width), names(component))
  width$extras <- sum(unlist(width)) - width$pchar
  width$term <- getOption("width")
  width$progress <- trunc((width$term - width$extras) / width$pchar)

  if (max <= min) stop("must have 'max' > 'min'")

  # Handles an update to the progress bar
  up <- function(value, extra = 0) {
    if (!is.finite(value) || value < min || value > max) {
      return()
    }
    .val <<- value
    nb <- round(width$progress * (value - min) / (max - min))
    pc <- round(100 * (value - min) / (max - min))
    extra <- round(100 * extra)
    if (nb == .nb && pc == .pc && .ex == extra) {
      return()
    }
    # Clear the current progress bar
    cat(paste0("\r", strrep(" ", width$term)))
    # Write the updated progress bar
    cat(paste0(
      "\r",
      component$prog_start,
      strrep(component$pchar, nb),
      strrep(" ", width$pchar * (width$progress - nb)),
      component$prog_end,
      sprintf(component$percent, pc),
      component$acc_sep,
      sprintf(component$acc_msg, extra)))
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
#' @param err_cond The original error condition that prompted this.
#'
#' @keywords internal
gibbs_step_err <- function(pmwgs, err_cond) {
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
  message("\nError while generating new group level parameters")
  message(err_cond)
  traceback()
  message("\nSaving current state of pmwgs object: ", sampler_tmp)
  saveRDS(pmwgs, file = sampler_tmp)
  stop("Stopping execution")
}
