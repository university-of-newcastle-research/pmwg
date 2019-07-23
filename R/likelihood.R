#' Function to calculate the log-likelihood of data, given random
#' effects, or coversely to produce synthetic data from given random
#' effects which match shape of "data". It relies on the rtdists package.
#' 
#' @param x A list of parameter values
#' @param data A data.frame containing variables for
#'        response time (rt), trial condition (condition)
#'        and accuracy (correct) for which the likelihood
#'        should be calculated.
#' @return The log likelihood of \code{x} given \code{data}.
#' @examples
#' x <- c(.2, .2, .2, .4, .3, 1.3, -2)
#' names(x) <- c("b1", "b2", "b3", "A", "v1", "v2", "t0")
#' rt <- c(0.432, 0.501, 0.31, 0.481, 0.342)
#' correct <- c(1, 2, 1, 1, 2)
#' condition <- c(1, 2, 3, 2, 1)
#' df <- data.frame(rt, correct, condition)
#' lba_loglike(x, df)
#' @export
lba_loglike <- function(x, data, sample = FALSE) {
  x <- exp(x)
  if (any(data$rt < x["t0"])) {
    return(-Inf)
  }
  # This is faster than "paste".
  bs <- x["A"] + x[c("b1", "b2", "b3")][data$condition]

  if (sample) {
    out <- rLBA(n = nrow(data), # nolint
                A = x["A"],
                b = bs,
                t0 = x["t0"],
                mean_v = x[c("v1", "v2")],
                sd_v = c(1, 1),
                distribution = "norm",
                silent = TRUE)
  } else {
    out <- dLBA(rt = data$rt, # nolint
                response = data$correct,
                A = x["A"],
                b = bs,
                t0 = x["t0"],
                mean_v = x[c("v1", "v2")],
                sd_v = c(1, 1),
                distribution = "norm",
                silent = TRUE)
    bad <- (out < 1e-10) | (!is.finite(out))
    out[bad] <- 1e-10
    out <- sum(log(out))
  }
  out
}
