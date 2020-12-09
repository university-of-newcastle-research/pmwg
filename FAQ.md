# Frequently Asked Questions (FAQ)

* [General Questions](#general-questions)
  * [I've found a problem with the package, what now?](#i've-found-a-problem-with-the-package,-what-now?)
* [Data Problems](#data-problems)
  * [My data are very complex, is there any way they can be stored as a list?](#my-data-are-very-complex,-is-there-any-way-they-can-be-stored-as-a-list?)
  * [I'm getting a very low acceptance rate, how can I check acceptance for each subject?](#i'm-getting-a-very-low-acceptance-rate,-how-can-i-check-acceptance-for-each-subject?)
## General Questions

### I've found a problem with the package, what now?

Head over to the [Issues](https://github.com/NewcastleCL/pmwg/issues) page and ask a question, suggest a new feature or report a bug.

## Data Problems

### My data are very complex, is there any way they can be stored as a list?

If your data is particularly complex and you currently have it saved as a list then the following conversion should move it into a `pmwg` compatible format.

```r
	# Play data - a list of data.frames in this case - but each list item should be able to be any object
	data_list <- list(tulving = data.frame(group = c('high', 'low'), meanrt = c(0.8, 1.4)), tversky = data.frame(group = c('high', 'low'), meanrt = c(0.9, 1.1)))
	# Turn that list into a tibble, which preserves the internal data.frame objects
	data <- tibble(subject = names(data_list), df = data_list)
	# Alternatively, use the I() function to interpret the list elements as is
	data <- data.frame(subject = names(data_list), df = I(data_list))
	# Then pass the data.frame/tibble as documented to the pmwgs object, and run your sampling stages.

	# The data passed into your log likelihood function would have the following form: (basically the output from either of the last two links above)
	# # A tibble: 1 x 2
	#   subject df              
	#   <chr>   <named list>    
	# 1 tulving <df[,2] [2 × 2]>
	# So you will also need to extract the internal data.frame or other object in your log likelihood function, using something like:
	df <- data$df[[data$subject]]
	# Then you can work with the data.frame (or whever format your list elements are) from there within the rest of your log-likelihood function.
```

### I'm getting a very low acceptance rate, how can I check acceptance for each subject?

If you are getting low acceptance rates (as displayed in the progress bar) sometimes this could be caused by one or two participant datasets.

The internal function `accept_rate` will return an array of mean acceptance rates, where the acceptance rate is the number of new, unique, particles accepted by the sampler

The default is to calculate this mean over the last 200 samples (or as many samples as exist if fewer than 200). The following two calls show how to do this for a custom windo, or for all samples.

```r
	# For a sampled object called sampler
	pmwg:::accept_rate(sampler, window_size=500)
	# For all samples
	pmwg:::accept_rate(sampler, window_size=sampler$samples$idx)
```
