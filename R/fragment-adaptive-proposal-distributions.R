
#'

## Test cases. Suppose these samples came from adaptation stage.
k <- 100 ## samples per random effect.

mu <- c(3, 3, 2, -1, 0)
sigma <- array(rnorm(length(mu)^2), dim = c(length(mu), length(mu)))
sigma <- sigma %*% t(sigma) ## Ensure posdef.
mu_samples <- t(rmvnorm(n = k, mean = mu, sigma = sigma))
sigma_samples <- array(sigma, dim = c(dim(sigma), k)) ## all "samples" identical.
sigma_samples <- sigma_samples * runif(length(sigma_samples), min = 0.99, max = 1.01)

## Turn the sigma samples into vectors via Cholesky + log
sigma_samples_unwound <- apply(sigma_samples, 3, unwind)

## These alpha samples are random effect draws, for just one participant.
alpha_samples <- rmvnorm(n = k, mean = mu, sigma = sigma)
alpha_samples <- array(t(alpha_samples), dim = c(length(mu), k))

## Put together the group and individual samples.
all_samples <- rbind(alpha_samples, mu_samples, sigma_samples_unwound)
## These are probably called mean_theta and covmat_theta in
## David's matlab code.
mu_tilde <- apply(all_samples, 1, mean)
sigma_tilde <- var(t(all_samples))

## Done. Here is how to use it to produce proposals.

## First, have a current proposal for group-level mu and alpha.
mu_proposal <- rmvnorm(1, mu, sigma)
sigma_proposal <- sigma_samples[, , 1] ## Just pick one for demo.

## These are the things that David calls "cond_var" and "cond_mean".
proposal_parameters <- condMVN(mean = mu_tilde, sigma = sigma_tilde, dependent.ind = 1:length(mu), given.ind = (length(mu) + 1):length(mu_tilde), X.given = c(mu_proposal, unwind(sigma_proposal)))

## Use them to sample as many proposals as you like.
proposals <- rmvnorm(n = 30, mean = proposal_parameters$condMean, sigma = proposal_parameters$condVar)
