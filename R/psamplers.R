#' Implementations of particle based sampling methods for model parameter estimation. Initially primarily an implementation of the Particle Metropolis within Gibbs sampler outlined in the paper available at https://arxiv.org/abs/1806.10089
#'
#' @docType package
#' @name psamplers
NULL

#
# Write functions only and document them with roxygen-styled comments.
# Example below taken from http://r-pkgs.had.co.nz/man.html
#

#' Add together two numbers.
#' 
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
#' @export
add <- function(x, y) {
      x + y
}
