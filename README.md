# psamplers #

<!-- badges: start -->
[![Build Status](https://travis-ci.org/newcastlecl/psamplers.svg?branch=master)](https://travis-ci.org/newcastlecl/psamplers)
[![Codecov test coverage](https://codecov.io/gh/NewcastleCL/samplers/branch/master/graph/badge.svg)](https://codecov.io/gh/NewcastleCL/samplers?branch=master)
<!-- badges: end -->


To install the currently recommended method is via devtools.

`devtools::install_github('newcastlecl/samplers')`

## Installing to an older version of R (< 3.6)

Terminal
* wget https://github.com/newcastlecl/samplers/archive/reduce_requirements.zip
* unzip reduce_requirements.zip
R
* devtools::install_version('mvtnorm', version='1.0-0')
* (Maybe need install.packages("mcmc")
* devtools::install_version('MCMCpack', version='1.4-0')
* install.packages(c("samplers-reduce_requirements/"), repos=NULL)
* 
