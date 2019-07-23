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
