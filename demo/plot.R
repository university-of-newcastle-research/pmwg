require(dplyr)
require(ggplot2)

plot.pmwgs <- function(x, type="mu", pars=NULL, subjects=NULL) {
  if (is.null(pars)) pars <- x$par_names
  if (is.null(subjects)) subjects <- x$subjects
  
  if (type == "mu") {
    data <- x$samples$theta_mu
    dimnames(data) <- list(x$par_names, NULL)
    data.frame(t(data)) %>%
      select(pars) %>%
      mutate(iter = 1:dim(data)[2]) %>%
      pivot_longer(-iter) %>%
      rename(parameter=name) %>%
      ggplot(mapping = aes(x=iter, y=value, col=parameter)) + geom_line()
  } else if (type == "alpha") {
    data <- x$samples$alpha
    dimnames(data) <- list("parameter" = x$par_names, "subject" = unique(x$data$subject), "iter" = 1:dim(data)[3])
    as.tbl_cube(data) %>%
      as_tibble() %>%
      dplyr::filter(parameter %in% pars) %>%
      dplyr::filter(subject %in% subjects) %>% 
      rename(value=data) %>% 
      ggplot(mapping = aes(x=iter, y=value, col=parameter)) + geom_line() + facet_wrap(~subject)
  } else if (type == "ll") {
    data <- x$samples$subj_ll
    dimnames(data) <- list("subject" = x$subjects, "iter" = 1:dim(data)[2])
    data.frame(t(data)) %>%
      select(subjects) %>%
      mutate(iter = 1:dim(data)[2]) %>%
      pivot_longer(-iter) %>%
      rename(parameter=name) %>%
      ggplot(mapping = aes(x=iter, y=value, col=parameter)) + geom_line()
  } else {
    stop("Unsupported plot type for pmwgs object")
  }
}
