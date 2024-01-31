# pmwg 0.2.7

## New features

* The `run_stage` function now takes a new argument: `p_accept`. The pmwg sampler object now stores a per subject epsilon (scale parameter for the covariance matrix). The `p_accept` argument is a target value for acceptance of new particles. Subjects with a low acceptance rate compared to the target lead to a lower epsilon and subjects with a high acceptance rate lead to a higher epsilon. Hyper parameters for this epsilon tuning come from: Garthwaite, P. H., Fan, Y., & Sisson, S. A. (2016).
* The `augment_sampler_epsilon` function allows the new pmwg package functions to work with older (saved) pmwg objects by extending the internal storage of the pmwg object and allowing subject specific epsilon values to be stored.
* Start points for sampling are now sampled from the prior rather than a standard normal. In cases where parameter values are not expected to be centred on zero this can lead to better start values.

## Minor improvements and fixes

* Improve error traceback when an error in the user supplied log-likelihood function is detected. Improves the ability of the enduser to debug their code.
* Improve progress bar formatting to remove warnings which caused display issues.
* BUGFIX: Calculation of scale parameter for a_half generation implemented incorrectly when sampling the mixture weights. Updated to reflect the original specification.
* Remove missing values from the debug object created when the gibbs step fails. Removes friction when debugging errors.
* Default on the number of particles for a `run_stage` call has been changed from 1000 to 100 particles. 100 particles should be more than sufficient for most models/data.
* Default on the number of particles for an `init` call has been changed from 1000 to 100 particles. 100 particles should be more than sufficient for most models/data.
* Default on the update frequency of the proposal distribution has been changed from 500 to 50 particles. This is a relatively cheap operation, so more frequent updates to improve estimation efficiency seems appropriate.
* Improve message formatting to make output cleaner and more informative.
* Better tests for malformed covariance matrixes have been imported into pmwg for cleaner reporting to the end user.
* Increases to the number of unique samples before attempting proposal distribution creation. With the new subject specific epsilon we are more likely to reach the `n_unique` value earlier. We now run for a bit longer to improve initial proposal distributions.
* Vendored in the random sampling from an inverse wishart distribution functions, allowing the removal of the MCMCpack dependency.
* Minor changes to documentation.

# pmwg 0.2.0

## New features

* The proposal distribution used in the sample stage will now be updated after a specified number of iterations, controlled by the new `run_stage` argument: `pdist_update_n`. By default this is set to 500, so after 500 iterations the proposal distribution will be regenerated and used for subsequent iterations of the sampler.
* The internal function `accept_rate` now operates over the entire sampler data store, and by default will return the rate of accepted newly generated random effects over a window covering the last 200 iterations. This will effect how the acceptance rate progress bar display progress, being more responsive to local changes in acceptance rate. This also means secondary runs of the adaptation stage can take into account unique random effect samples from previous runs for it's internal testing.
* New function, `relabel_samples`, added to allow a user to change the label indicating the samples stage of origin. This function takes a pmwg samplers object, a set of indices corresponding to the samples you wish to relabel and two arguments detailing to expected value for the current stage (default "burn") and the stage to would like to relabel the samples to (default "adapt"). The goal for the function is to take samples generated in a burn in stage that after visualising appear to have successfully converged around the posterior and make them available for internal function that check for adaptation and generate the proposal distribution used in the final sampling stage.
* New included data object `sampled_forstmann`, a pmwgs object with a simple LBA model (requires rtdists package) applied to the previously included `forstmann` dataset. Includes 3 runs of the sampler with low numbers of iterations and particles. It has the data element stripped so see the Details in the documentation for the object (`?sampled_forstmann`) for advice on using the object.

## Minor improvements and fixes

* BUGFIX: If all random effects were unique the previous version would not test the ability to create the proposal distribution and leave adaptation stage early. Now fixed by explicitly splitting the first parameter matrix to a list before testing.
* Improved error messages from the `gibbs_step` internal function
* Changed Acc to New in progress bar to focus that it is the rate of newly generated particles being accepted.
* Documentation improvements:
  * Addition of some examples using the sampled_forstmann object
  * Updates to accept_progress_bar documentation
  * Minor edits to documentation elsewhere.
* New dependency `checkmate` used for testing arguments to `run_stage`

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
