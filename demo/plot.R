require(dplyr)
require(ggplot2)
require(tidyr)

plot.pmwgs <- function(x, type="mu", pars=NULL, subjects=NULL, iters=NULL, transform.func=NULL) {
  if (is.null(pars)) pars <- x$par_names
  if (is.null(subjects)) subjects <- x$subjects
  if (is.null(iters)) iters <- 1:x$samples$idx

  if (type == "mu") {
    data <- x$samples$theta_mu[, iters]
    if (!is.null(transform.func)) data <- transform.func(data)
    dimnames(data) <- list(x$par_names, NULL)
    data.frame(t(data)) %>%
      select(pars) %>%
      mutate(iter = 1:dim(data)[2]) %>%
      pivot_longer(-iter) %>%
      rename(parameter = name) %>%
      ggplot(mapping = aes(x = iter, y = value, col = parameter)) + geom_line()
  } else if (type == "alpha") {
    data <- x$samples$alpha[, , iters]
    if (!is.null(transform.func)) data <- transform.func(data)
    dimnames(data) <- list("parameter" = x$par_names,
                           "subject" = unique(x$data$subject),
                           "iter" = 1:dim(data)[3])
    as.tbl_cube(data) %>%
      as_tibble() %>%
      dplyr::filter(parameter %in% pars) %>%
      dplyr::filter(subject %in% subjects) %>%
      rename(value = data) %>%
      ggplot(mapping = aes(x = iter, y = value, col = parameter)) +
        geom_line() +
        facet_wrap(~subject)
  } else if (type == "ll") {
    data <- x$samples$subj_ll[, iters]
    if (!is.null(transform.func)) data <- transform.func(data)
    data <- t(data)
    dimnames(data) <- list("iter" = 1:dim(data)[1], "subject" = x$subjects)
    data.frame(data, check.names=FALSE) %>%
      mutate(iter = 1:dim(data)[1]) %>%
      pivot_longer(-iter) %>%
      rename(subject=name) %>%
      dplyr::filter(subject %in% subjects) %>%
      ggplot(mapping = aes(x = iter, y = value, col = subject)) + geom_line()
  } else {
    stop("Unsupported plot type for pmwgs object")
  }
}
