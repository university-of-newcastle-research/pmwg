# psamplers #

[![Build Status](https://travis-ci.org/gjcooper/psamplers.svg?branch=master)](https://travis-ci.org/gjcooper/psamplers)
[![codecov.io](https://codecov.io/github/gjcooper/psamplers/coverage.svg?branch=master)](https://codecov.io/github/gjcooper/psamplers?branch=master)

To install the currently recommended method is via devtools.

`devtools::install_github('gjcooper/samplers')`

## Installing to an older version of R (< 3.6)

Terminal
* wget https://github.com/gjcooper/samplers/archive/reduce_requirements.zip
* unzip reduce_requirements.zip
R
* devtools::install_version('mvtnorm', version='1.0-0')
* (Maybe need install.packages("mcmc")
* devtools::install_version('MCMCpack', version='1.4-0')
* install.packages(c("samplers-reduce_requirements/"), repos=NULL)
* 
