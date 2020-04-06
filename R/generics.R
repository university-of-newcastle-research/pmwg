#' Initialise a sampler with starting points
#'
#' Takes the starting values for the individual and possibly group
#' level parameters, in preparation for further iterations,
#' whether they be burn-in or sampling.
#'
#' @param x The sampling object that provides the parameters.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The sampler object but with sampled values for subject and group
#'   parameters. Exact sampling method depends on the stage being run.
#' @examples
#' # No example yet
#' @export
run_stage <- function(x, ...) {
  if (is.null(attr(x, "class"))) {
    print("No object to run a stage on")
  }
  else UseMethod("run_stage")
}


#' Initialise a sampler with starting points
#'
#' Takes the starting values for the individual and possibly group
#' level parameters, in preparation for further iterations,
#' whether they be burn-in or sampling.
#'
#' @param x The sampling object that provides the parameters.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The sampler object but with initial values set for all parameter
#'   types
#' @examples
#' # No example yet
#' @export
init <- function(x, ...) {
  if (is.null(attr(x, "class"))) {
    print("No object to add start points to")
  }
  else UseMethod("init")
}
