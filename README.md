# pmwg - Particle Metropolis within Gibbs <img src="man/figures/hexlogo_small.png" align="right"/> #

<!-- badges: start -->
[![R build status](https://github.com/NewcastleCL/pmwg/workflows/R-CMD-check/badge.svg)](https://github.com/NewcastleCL/pmwg/actions)
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
The document there includes the motivation for the approach, several detailed examples of the package in action and a list of common problems and troubleshooting techniques.
Also available online is the package documentation at https://newcastlecl.github.io/pmwg/ which consists of this README, a Reference of help documentation for individual functions, a list of changes to the project over time and more.
Finally there is a page containing some frequently asked questions which can be found at https://newcastlecl.github.io/pmwg/FAQ.html

Included on the pmwg website is also a getting started guide to the package, available from https://newcastlecl.github.io/pmwg/articles/pmwg.html

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
