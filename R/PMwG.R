# Function to calculate the log-likelihood of data, given random
# effects, or coversely to produce synthetic data from given random
# effects which match shape of "data".
ll <- function(x, data, sample = FALSE) {
  x <- exp(x)
  if (any(data$rt < x["t0"])) {
    return(-Inf)
  }
  # b.ind=paste0("b",data$condition)
  # bs=x["A"]+x[b.ind]
  bs <- x["A"] + x[c("b1", "b2", "b3")][data$condition] # This is faster than "paste".
  if (sample) {
    out <- rLBA(n = nrow(data), A = x["A"], b = bs, t0 = x["t0"], mean_v = x[c("v1", "v2")], sd_v = c(1, 1), distribution = "norm", silent = TRUE)
  } else {
    out <- dLBA(rt = data$rt, response = data$correct, A = x["A"], b = bs, t0 = x["t0"], mean_v = x[c("v1", "v2")], sd_v = c(1, 1), distribution = "norm", silent = TRUE)
    bad <- (out < 1e-10) | (!is.finite(out))
    out[bad] <- 1e-10
    out <- sum(log(out))
  }
  out
}

# Paralellisable particle sampling function (cond/orig/real/log).
pm.corl <- function(s, data, n.particles, mu, sig2, particles) {
  # This uses the simplest, and slowest, proposals: mixture of the
  # the popultion distribution and gaussian around current random effect.
  wmix <- 0.5
  n1 <- rbinom(n = 1, size = n.particles, prob = wmix)
  if (n1 < 2) n1 <- 2
  if (n1 > (n.particles - 2)) n1 <- n.particles - 2 # These just avoid degenerate arrays.
  proposals1 <- rmvnorm(n1, mu, sig2)
  proposals2 <- rmvnorm(n.particles - n1, particles[, s], sig2)
  proposals <- rbind(proposals1, proposals2)
  colnames(proposals) <- names(mu) # stripped otherwise.
  proposals[1, ] <- particles[, s] # Put the current particle in slot 1.
  lw <- apply(proposals, 1, ll, data = data[data$subject == s, ]) # Density of data given random effects proposal.
  lp <- dmvnorm(x = proposals, mean = mu, sigma = sig2, log = TRUE) # Density of random effects proposal given population-level distribution.
  lm <- log(wmix * exp(lp) + (1 - wmix) * dmvnorm(x = proposals, mean = particles[, s], sigma = sig2)) # Density of proposals given proposal distribution.
  l <- lw + lp - lm # log of importance weights.
  weight <- exp(l - max(l))
  proposals[sample(x = n.particles, size = 1, prob = weight), ]
}
