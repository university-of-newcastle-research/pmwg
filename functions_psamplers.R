
## Creat pp.data

generate.ppdata <- function (sampled,n.posterior=NULL,stage="sample",ll_func = sampled$ll_func) {
  index = which(sampled$samples$stage==stage)
  if (is.null(n.posterior)){
    iterations=round(seq(from=index[1],to=index[length(index)],length.out=length(index)))
  }
  else {
    iterations=round(seq(from=index[1],to=index[length(index)],length.out=n.posterior))
  }
  data = split(x=sampled$data,f=sampled$data$subject)
  S = unique(sampled$data$subject)
  pp.data=list()
  for (s in S){
    for (i in 1:length(iterations)) {
      x <- sampled$samples$alpha[,s,iterations[i]]
      names(x) <- pars
      out <- ll_func(x=x,data=data[[s]],sample=TRUE)
      if (i==1){
        pp.data[[s]] = cbind(pp.subject=i,out)
      }
      else {
        pp.data[[s]]= rbind(pp.data[[s]],cbind(pp.subject=i,out))
      }
    }
  }
  names(pp.data)=names(data)
  tidy.pp.data=do.call(rbind,pp.data)
  return(tidy.pp.data)
}


hist_accept <- function(x, stage, ...) {
  main.name <- switch(stage,
                      burn = "Burnin",
                      adapt = "Adaptation",
                      sample = "Sampled")
  hist(x,
       main = paste("Acceptance Rate", main.name),
       xlab = "Acceptance Rate",
       ...)
  abline(v = mean(x), col = "red")
}

acceptance_rate = function(x, stage = NULL,plot = FALSE, ...) {
  stages = unique(x$samples$stage)
  no.init = "init"
  stages = stages [!stages %in% no.init]
  tmp1 <-
    lapply(setNames(stages, stages), function(i)
      x$samples$alpha[, , x$samples$stage == i])
  tmp2 = lapply(tmp1, function(i)
    apply(i[1, , ], 1, diff))
  tmp3 = lapply(tmp2, function(i)
    apply(i != 0, 2, mean))
  
  if (!is.null(stage)) {
    tmp35 <- tmp3[stage]
    if (plot) {
      op <- par (mfrow = c(1, 1))
      tmp4 = unname(unlist(tmp35))
      hist_accept(tmp4, stage, ...)
      par(op)
    }
    class(tmp35) <- "acc_rate"
    return(tmp35)
  }
  else {
    if (plot) {
      op <- par (mfrow = c(1, length(tmp3)))
      lapply(seq_along(tmp3), function(i)
        hist_accept(tmp3[[i]], stage = stages[i], ...))
      par(op)
    }
    class(tmp3) <- "acc_rate"
    return(tmp3)
  }
  
}


#calculate summary stats
summary.acc_rate <- function(object, what, ...) {
  # what needs to be a character vector with function names that can be applied
  # to a numeric vector (e.g., "mean", "min", "max", etc.)
  out <- lapply(setNames(what, what), function(i) lapply(object, i))
  class(out) <- "summary.acc_rate"
  return(out)
}
