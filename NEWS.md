# pmwg 0.1.9

* Clean up documentation ready for CRAN submission
* Add extra argument to as_mcmc to filter specific samples

# pmwg 0.1.8

* Fixes: Non-zero mean in prior was not applied in the sampling algorithm
* Updated author/creator information.

# pmwg 0.1.7

* Inform users of the default epsilon value chosen by the package if no explicit value is given.
* Rename list elements for the prior distribution on theta_mu to increase clarity over their purpose.

# pmwg 0.1.6

* Add `as_mcmc` function to export mcmc objects (from coda package) or a list of mcmc objects for the covariance matrix and subject random effects.
* Add unique subject identifiers as the appropriate dimnension names in sample storage (ie for random effects)
* Fixes: dimension names no longer lost when updating the sample storage with new samples.
* New particles argument in the `init` function to reduce time for sampler initialisation.
 
# pmwg 0.1.5

* Closes Issue #16 - Better errors on incorrect arguments to run_stage
* Closes Issue #26 - Detect and raise error on missing subject column in data frame 

# pmwg 0.1.4

* Update package name internally to pmwg
* Separate adaptation check code to new function

# pmwg 0.1.2

* Added a `NEWS.md` file to track changes to the package.
* Initial version in preparation for CRAN submission
