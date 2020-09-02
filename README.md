# pmwg - Particle Metropolis within Gibbs <img src="man/figures/hexlogo_small.png" align="right"/> #

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/NewcastleCL/pmwg.svg?branch=release)](https://travis-ci.com/NewcastleCL/pmwg)
[![Codecov test coverage](https://codecov.io/gh/NewcastleCL/pmwg/branch/release/graph/badge.svg)](https://codecov.io/gh/NewcastleCL/pmwg?branch=release)
<!-- badges: end -->

## Installation

To install the latest stable version you can use the following command to install from CRAN:

`install.packages("pmwg")`

If you want the (possibly unstable) development version, you can also install the package using devtools as follows:

`devtools::install_github('newcastlecl/pmwg', ref="develop")`

This package is tested and should work on all versions of R > 3.5, however instructions on installing to an earlier version of R are included below.

## Using the package

The best introduction to the package can be found at the bookdown site located at: https://newcastlecl.github.io/samplerDoc/
The documentation there includes the motivation for the approach, several detailed examples of the package in action and a list of common problems and troubleshooting techniques. Included here is a skeleton of the required steps to run the sampler.

### Quickstart Guide

```r
library(pmwg)

# Create a log likelikhood function for your model
loglikelihood_func <- function(x, data) {
  # x are the current model parameter values
  # return the log likelihood for the data (for the subject) given x
}

# Create the sampler object with your data, parameter names, loglike function and
# A list of priors for theta_mu_mean (mean of model parameters) and theta_mu_var (covariance of model parameters)
sampler <- pmwgs(
  data = my_data_source,
  pars = c("a", "list", "of", "par", "names"),
  ll_func = loglikelihood_func,
  prior = list(theta_mu_mean = rep(0, length(pars)), theta_mu_var = diag(rep(1, length(pars))))
)

# Initialise (generate first random effects for sampler)
sampler <- init(sampler)  # Can also pass start points for sampler

# Run each stage of the sampler, can adjust number of particles on each
sampler <- run_stage(sampler, stage="burn")
sampler <- run_stage(sampler, stage="adapt")
sampler <- run_stage(sampler, stage="sample")
```

The `run_stage` command can also be passed other arguments such as `iter` for number of iterations, `particles` for number of particles among others. For a full list see [the description in the PMwG Tutorial Book](https://newcastlecl.github.io/samplerDoc/pmwg-sampler-and-signal-detection-theory.html#run-sdtsampler).

## Installing to an older version of R (< 3.5)

**Terminal**

```bash
wget https://github.com/newcastlecl/pmwg/archive/reduce_requirements.zip
unzip reduce_requirements.zip
```

**R**

```R
devtools::install_version('mvtnorm', version='1.0-0')
# (Maybe need install.packages("mcmc")
devtools::install_version('MCMCpack', version='1.4-0')
install.packages(c("pmwg-reduce_requirements/"), repos=NULL)
```
