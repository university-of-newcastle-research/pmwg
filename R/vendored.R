##
## Wishart and Inverse wishart distribution functions vendored from MCMCpack
## version 1.6-3
##
## Wishart based on code originally posted by Bill Venables to S-news
## on 6/11/1998
##
## KQ on 2/5/2001
##
## GC added to pmwg 15/01/2024

# rwish delivers a pseudo-random Wishart deviate
#
# USAGE:
#
#   A <- rwish(v, S)
#
# INPUT:
#
#   v    degrees of freedom
#
#   S    Scale matrix
#
# OUTPUT:
#
#  A     a pseudo-random Wishart deviate
#
# riwish generates a draw from the inverse Wishart distribution
# (using the Wishart generator)

#' The Wishart Distribution
#'
#' Random generation from the Wishart distribution.
#'
#' The mean of a Wishart random variable with \code{v} degrees of freedom and
#' inverse scale matrix \code{S} is \eqn{vS}.
#'
#' @aliases dwish rwish
#'
#' @param v Degrees of freedom (scalar).
#'
#' @param S Inverse scale matrix \eqn{(p \times p)}.
#'
#' @return \code{rwish} generates one random draw from the distribution.
#'
#' @keywords internal
#'
#' @examples \dontrun{
#' draw <- rwish(3, matrix(c(1,.3,.3,1),2,2))
#' }
rwish <- function(v, S) {
	if (!is.matrix(S))
		S <- matrix(S)
	if (nrow(S) != ncol(S)) {
		stop(message="S not square in rwish().\n")
	}
	if (v < nrow(S)) {
		stop(message="v is less than the dimension of S in rwish().\n")
	}
	p <- nrow(S)
	CC <- chol(S)
	Z <- matrix(0, p, p)
	diag(Z) <- sqrt(stats::rchisq(p, v:(v-p+1)))
	if(p > 1) {
		pseq <- 1:(p-1)
		Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- stats::rnorm(p*(p-1)/2)
	}
	return(crossprod(Z %*% CC))
}

##
## Inverse Wishart
##

#' The Inverse Wishart Distribution
#'
#' Random generation from the Inverse Wishart distribution.
#'
#' The mean of an inverse Wishart random variable with \code{v} degrees of
#' freedom and scale matrix \code{S} is \eqn{(v-p-1)^{-1}S}.
#'
#' @param v Degrees of freedom (scalar).
#'
#' @param S Scale matrix \eqn{(p \times p)}.
#'
#' @return \code{riwish} generates one random draw from the distribution.
#'
#' @keywords internal
#'
#' @examples \dontrun{
#' draw <- riwish(3, matrix(c(1,.3,.3,1),2,2))
#' }
riwish <- function(v, S) {
	return(solve(rwish(v,solve(S))))
}
