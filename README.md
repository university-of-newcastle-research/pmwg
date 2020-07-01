# pmwg - Particle Metropolis within Gibbs Sampler #

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/NewcastleCL/pmwg.svg?branch=release)](https://travis-ci.com/NewcastleCL/pmwg)
[![Codecov test coverage](https://codecov.io/gh/NewcastleCL/pmwg/branch/release/graph/badge.svg)](https://codecov.io/gh/NewcastleCL/pmwg?branch=release)
<!-- badges: end -->

## Installation

To install the currently recommended method is via devtools.

`devtools::install_github('newcastlecl/pmwg', ref="release")`

## Using the package

### Quickstart Guide

```r
library(pmwg)

# Create a log likelikhood function for your model
loglikelihood_func <- function(x, data) {
  # x are the current model parameter values
  # return the log likelihood for the data (for the subject) given x
}

# Create the sampler object with your data, parameter names, loglike function and
# A list of priors for theta_mu (model parameters) and theta_sig (covariance matrix)
sampler <- pmwgs(
  data = my_data_source,
  pars = c("a", "list", "of", "par", "names"),
  ll_func = loglikelihood_func,
  prior = list(theta_mu = rep(0, length(pars)), theta_sig = diag(rep(1, length(pars))))
)

# Initialise (generate first random effects for sampler)
sampler <- init(sampler)  # Can also pass start points for sampler

# Run each stage of the sampler, can adjust number of particles on each
sampler <- run_stage(sampler, stage="burn")
sampler <- run_stage(sampler, stage="adapt")
sampler <- run_stage(sampler, stage="sample")
```

The `run_stage` command can also be passed other arguments such as `iter` for number of iterations, `particles` for number of particles among others. For a full list see: https://newcastlecl.github.io/samplerDoc/pmwg-sampler-and-signal-detection-theory.html#run-sdtsampler


The documentation is available at https://newcastlecl.github.io/samplerDoc/.


## Installing to an older version of R (< 3.6)

Terminal
* wget https://github.com/newcastlecl/pmwg/archive/reduce_requirements.zip
* unzip reduce_requirements.zip
R
* devtools::install_version('mvtnorm', version='1.0-0')
* (Maybe need install.packages("mcmc")
* devtools::install_version('MCMCpack', version='1.4-0')
* install.packages(c("pmwg-reduce_requirements/"), repos=NULL)
