  do.call("rbind",events),
  purrr::map_df(events,dplyr::bind_rows),
  times=1000)
ggplot2::autoplot(mb)+theme_bw()
ggplot2::autoplot(mb)#+theme_bw()
mb <- microbenchmark(
#  plyr::rbind.fill(events),
  dplyr::bind_rows(events),
  data.table::rbindlist(events),
  plyr::ldply(events,data.frame),
  do.call("rbind",events),
  purrr::map_df(events,dplyr::bind_rows),
  Reduce(rbind, events)
  times=1000)
mb <- microbenchmark(
#  plyr::rbind.fill(events),
  dplyr::bind_rows(events),
  data.table::rbindlist(events),
  plyr::ldply(events,data.frame),
  do.call("rbind",events),
  purrr::map_df(events,dplyr::bind_rows),
  Reduce(rbind, events),
  times=1000)
mb <- microbenchmark(
  plyr::rbind.fill(events),
#  dplyr::bind_rows(events),
  data.table::rbindlist(events),
  plyr::ldply(events,data.frame),
  do.call("rbind",events),
  purrr::map_df(events,dplyr::bind_rows),
  Reduce(rbind, events),
  times=1000)
?dplyr::bind_rows
?data.table::rbindlist
mb <- microbenchmark(
  plyr::rbind.fill(events),
  dplyr::bind_rows(events),
  data.table::rbindlist(events, use.names=TRUE),
  plyr::ldply(events,data.frame),
  do.call("rbind",events),
  purrr::map_df(events,dplyr::bind_rows),
  Reduce(rbind, events),
  times=1000)
ggplot2::autoplot(mb)#+theme_bw()
?identical
events <- lapply(raw_data, FUN = function(data_line) {
  df <- read.csv(text = data_line$jsPsychData)
})
evtest0 <- plyr::rbind.fill(events)
evtest1 <- dplyr::bind_rows(events)
evtest2 <- data.table::rbindlist(events, use.names=TRUE)
evtest3 <- plyr::ldply(events,data.frame)
evtest4 <- do.call("rbind",events)
evtest5 <- purrr::map_df(events,dplyr::bind_rows)
evtest6 <- Reduce(rbind, events)
?dplyr::all_equal
dplyr::all_equal(evtest0, evtest1)
dplyr::all_equal(evtest0, evtest2)
dplyr::all_equal(evtest0, evtest3)
dplyr::all_equal(evtest0, evtest4)
dplyr::all_equal(evtest0, evtest5)
dplyr::all_equal(evtest0, evtest6)
??dplyr::bind_rows
?dplyr::bind_rows
warnings()
gplot2::autoplot(mb)#+theme_bw()
str(events)
warnings()
dplyr::all_equal(evtest0, evtest6)
dplyr::all_equal(evtest0, evtest5)
str(events)
events <- lapply(raw_data, FUN = function(data_line) {
  df <- read.csv(text = data_line$jsPsychData)
})
events <- Reduce(rbind, events)
str(events)
as.character(events$rt) == evtest1$rt
all(as.character(events$rt) == evtest1$rt)
NvimR.paragraph()
NvimR.paragraph()
pref <- events[events$type == "preference", pref_columns]
pref$rt <- as.numeric(as.character(pref$rt))
mds <- events[events$type == "rating", mds_columns]
mds$rt <- as.numeric(as.character(mds$rt))
head(pref)
str(pref)
unique(pref$attr_order)
pref[pref$pair_id == 4839, ]
NvimR.paragraph()
library(ggplot2)
NvimR.paragraph()
NvimR.paragraph()
NvimR.paragraph()
NvimR.paragraph()
NvimR.paragraph()
str(pref)
pref$sonaID <- as.factor(pref$sonaID)
NvimR.paragraph()
pref$rt <- as.numeric(as.character(pref$rt))/1000
pref$sonaID <- as.factor(pref$sonaID)
mds <- events[events$type == "rating", mds_columns]
mds$rt <- as.numeric(as.character(mds$rt))
write.csv(pref, file = "pref_data.csv")
write.csv(mds, file = "mds_data.csv")
ggplot(pref, aes(x=rt, color=sonaID)) +
  geom_histogram(fill="white", alpha=0.5, position="identity")
library(patchwork)
NvimR.paragraph()
library(dplyr)
pref %>% group_by(sonaID) %>% summarize()
pref %>% group_by(sonaID) %>% summarize(grp.mean=mean(rt))
ibrary(dplyr)
library(dplyr)
NvimR.paragraph()
NvimR.paragraph()
NvimR.paragraph()
NvimR.paragraph()
pref_hist / mds_hist
NvimR.paragraph()
NvimR.paragraph()
NvimR.paragraph()
NvimR.paragraph()
NvimR.paragraph()
NvimR.paragraph()
NvimR.paragraph()
NvimR.paragraph()
pref_hist / mds_hist
NvimR.paragraph()
NvimR.paragraph()
pref_hist / mds_hist
NvimR.paragraph()
NvimR.paragraph()
NvimR.paragraph()
pref_hist / mds_hist
NvimR.paragraph()
NvimR.paragraph()
pref_hist / mds_hist
NvimR.paragraph()
NvimR.paragraph()
pref_hist / mds_hist
NvimR.paragraph()
NvimR.paragraph()
pref_hist / mds_hist
quit(save = "no")
?png
:vsp
vsp
devtools::install_github("gaborcsardi/prompt")
?sample
source("~/.Rprofile")
source("~/.Rprofile")
ls
ls()
devtools::check()
usethis::use_testthat()
diag(rep(1,7))
rnorm(4*4, dim=c(4,4))
rnorm(4*4)
?rnorm
?matrix
matrix(rnorm(16), nrow=4, ncol=4)
devtools::test()
devtools::test()
devtools::test()
devtools::test()
test_mat <- matrix(rnorm(16), nrow=4, ncol=4)
test_mat2 <- test_mat %*% t(test_mat)
test_mat
test_mat2
psampelrs::unwind(test_mat)
psamplers::unwind(test_mat)
psamplers::unwind(test_mat2)
lower.tri(test_mat2)
test_mat2[lower.tri(test_mat2)]
test_mat3 <- diag(rep(1, 7))
psamplers::unwind(test_mat3)
test_mat3 <- diag(rep(1, 4))
psamplers::unwind(test_mat3)
psamplers::wind(psamplers::unwind(test_mat3))
psamplers::wind(rnorm(10))
psamplers::unwind(psamplers::wind(rnorm(10)))
dim(psamplers::unwind(psamplers::wind(rnorm(10))))
dim(psamplers::wind(rnorm(10)))
dim(psamplers::unwind(psamplers::wind(rnorm(10))))
length(psamplers::unwind(psamplers::wind(rnorm(10))))
type(length(psamplers::unwind(psamplers::wind(rnorm(10)))))
class(length(psamplers::unwind(psamplers::wind(rnorm(10)))))
devtools::test()
devtools::test()
?expect_error
devtools::test()
devtools::test()
usethis::use_coverage()
usethis::use_coverage()
usethis::use_coverage
usethis::use_coverage()
usethis::use_coverage()
ls
usethis::use_coverage()
usethis::use_coverage()
ls
:qa
usethis::use_coverage()
usethis::coverage
usethis::use_coverage
usethis::use_codecov_badge
usethis::use_coveralls_badge
?glue
usethis::use_coverage()
usethis::use_travis()
ls()
load("~/Downloads/wagenmakers2008.RData")
ls()
str(data)
library(rtdists)
str(speed_acc)
speed_acc
speed_acc
head(data)
wgnmks2008Fast <- as.data.frame(table(wgnmks2008$subject, wgnmks2008$cond,
                                  wgnmks2008$stim, wgnmks2008$resp))
names(wgnmks2008Fast) <- c("subject", "cond", "stim", "resp", "n")
wgnmks2008 <- data
wgnmks2008Fast <- as.data.frame(table(wgnmks2008$subject, wgnmks2008$cond,
                                  wgnmks2008$stim, wgnmks2008$resp))
head(data)
head(data,10)
str(data)
wgnmks2008 <- data.frame(subject = data$subject, cond = data$prop)
str(wgnmks2008)
wgnmks2008 <- data.frame(subject = data$subject, cond = data$prop, stim = data$freq, resp = data$resp, rt = data$rt, correct = data$correct))
wgnmks2008 <- data.frame(subject = data$subject, cond = data$prop, stim = data$freq, resp = data$resp, rt = data$rt, correct = data$correct)
head(wgnmks2008)
SDT_ll <- function(x, data, sample=FALSE){
  if (sample){
    data$response <- NA
  } else {
    out <- numeric(nrow(data))
  }
  if (!sample){
  for (i in 1:nrow(data)) {
    if (data$cond[i] == "w") {
    if (data$stim[i] == "hf") {
      if (data$resp[i] == "W") {
        out[i] <- pnorm(x["C.w"], mean = x["HF.d"], sd = 1,
                       log.p = TRUE, lower.tail = FALSE)
      } else {
        out[i] <- pnorm(x["C.w"], mean = x["HF.d"], sd = 1,
                       log.p = TRUE, lower.tail = TRUE)
      }
    } else if (data$stim[i] == "lf"){
      if (data$resp[i] == "W"){
        out[i] <- pnorm(x["C.w"], mean = x["LF.d"], sd = 1,
                       log.p = TRUE, lower.tail = FALSE)
      } else {
        out[i] <- pnorm(x["C.w"], mean = x["LF.d"], sd = 1,
                       log.p = TRUE, lower.tail = TRUE)
        }
      } else if (data$stim[i] == "vlf") {
        if (data$resp[i] == "W") {
          out[i] <- pnorm(x["C.w"], mean = x["VLF.d"], sd = 1,
                         log.p = TRUE, lower.tail = FALSE)
        } else {
          out[i] <- pnorm(x["C.w"], mean = x["VLF.d"], sd = 1,
                         log.p = TRUE, lower.tail = TRUE)
          }
        } else {
      if (data$resp[i] == "W") {
        out[i] <- pnorm(x["C.w"], mean = 0, sd = 1,
                       log.p = TRUE, lower.tail = FALSE)
      } else {
        out[i] <- pnorm(x["C.w"], mean = 0, sd = 1,
                       log.p = TRUE, lower.tail = TRUE)
        }
      }
    } else {
      if (data$stim[i] == "hf") {
        if (data$resp[i] == "W") {
          out[i] <- pnorm(x["C.nw"], mean = x["HF.d"], sd = 1, 
                         log.p = TRUE, lower.tail = FALSE)
        } else {
          out[i] <- pnorm(x["C.nw"], mean = x["HF.d"], sd = 1, 
                         log.p = TRUE, lower.tail = TRUE)
        }
      }  else if (data$stim[i] == "lf") {
        if (data$resp[i] == "W") {
          out[i] <- pnorm(x["C.nw"], mean = x["LF.d"], sd = 1, 
                         log.p = TRUE, lower.tail = FALSE)
        } else {
          out[i] <- pnorm(x["C.nw"], mean = x["LF.d"], sd = 1, 
                         log.p = TRUE, lower.tail = TRUE)
        }
      } else if (data$stim[i] == "vlf") {
        if (data$resp[i] == "W") {
          out[i] <- pnorm(x["C.nw"], mean = x["VLF.d"], sd = 1, 
                         log.p = TRUE, lower.tail = FALSE)
        } else {
          out[i] <- pnorm(x["C.nw"], mean = x["VLF.d"], sd = 1, 
                         log.p = TRUE, lower.tail = TRUE)
        }
      } else {
        if (data$resp[i] == "W") {
          out[i] <- pnorm(x["C.nw"], mean = 0, sd = 1, 
                         log.p = TRUE, lower.tail = FALSE)
        } else {
          out[i] <- pnorm(x["C.nw"], mean = 0, sd = 1, 
                         log.p = TRUE, lower.tail = TRUE)
        }
      }
    }
  }
  sum(out)
  }
}
ls()
pars <- log(c(C.w = 1, C.nw = 0.5, HF.d = 3, LF.d = 1.8, VLF.d = 0.7))
SDT_ll(pars, wgnmks2008, sample = FALSE)
wgnmks2008Fast <- as.data.frame(table(wgnmks2008$subject, wgnmks2008$cond,
                                  wgnmks2008$stim, wgnmks2008$resp))
names(wgnmks2008Fast) <- c("subject", "cond", "stim", "resp", "n")
SDT_ll_fast <- function(x, data, sample = FALSE) {
  if (!sample) {
    out <- numeric(nrow(data))
    for (i in 1:nrow(data)) {
      if (data$cond[i] == "w") {
        if (data$stim[i] == "hf") {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["HF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["HF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
            }
          } else if (data$stim[i] == "lf") {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["LF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["LF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        } else if (data$stim[i] == "vlf") {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["VLF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["VLF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        } else {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = 0,
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = 0,
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        }
      }else{
        if (data$stim[i] == "hf") {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["HF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["HF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        }  else if (data$stim[i] == "lf"){
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["LF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["LF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        } else if (data$stim[i] == "vlf") {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["VLF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["VLF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        } else {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = 0,
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = 0,
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        }
      }
    }
    sum(out)
  }
ls()
SDT_ll_fast <- function(x, data, sample = FALSE) {
  if (!sample) {
    out <- numeric(nrow(data))
    for (i in 1:nrow(data)) {
      if (data$cond[i] == "w") {
        if (data$stim[i] == "hf") {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["HF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["HF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
            }
          } else if (data$stim[i] == "lf") {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["LF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["LF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        } else if (data$stim[i] == "vlf") {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["VLF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["VLF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        } else {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = 0,
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.w"], mean = 0,
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        }
      }else{
        if (data$stim[i] == "hf") {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["HF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["HF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        }  else if (data$stim[i] == "lf"){
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["LF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["LF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        } else if (data$stim[i] == "vlf") {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["VLF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["VLF.d"],
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        } else {
          if (data$resp[i] == "W") {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = 0,
                                     sd = 1, log.p = TRUE, lower.tail = FALSE)
          } else {
            out[i] <- data$n[i] * pnorm(x["C.nw"], mean = 0,
                                     sd = 1, log.p = TRUE, lower.tail = TRUE)
          }
        }
      }
    }
    sum(out)
  }
}
ls()
SDT_ll(pars, wgnmks2008, sample = FALSE)
SDT_ll_fast(pars, wgnmks2008Fast, sample = FALSE)
pars <- c("C.w","C.nw","HF.d","LF.d","VLF.d") # This is the same as the `pars` vector specified above
priors <- list(
  theta_mu = rep(0, length(pars)),
  theta_sig = diag(rep(1, length(pars)))
)
sampler <- pmwgs(
  data = wgnmks2008Fast,
  pars = pars,
  prior = priors,
  ll_func = SDT_ll_fast
)
devtools::load_all()
sampler <- pmwgs(
  data = wgnmks2008Fast,
  pars = pars,
  prior = priors,
  ll_func = SDT_ll_fast
)
start_points <- list(
  mu = rep(0, length.out = length(pars)),
  sig2 = diag(rep(.01, length(pars)))
)
sampler <- init(sampler, theta_mu = start_points$mu,
                theta_sig = start_points$sig2)
burned <- run_stage(sampler, stage = "burn", iter = 1000, particles = 1000, display_progress = TRUE, n_cores = 4); adapted <- run_stage(burned, stage = "adapt", iter = 1000, particles = 1000, n_cores = 4); sampled <- run_stage(adapted, stage = "sample", iter = 1000, particles = 100, n_cores = 4)
burned <- run_stage(sampler, stage = "burn", iter = 300, particles = 100, display_progress = TRUE, n_cores = 4); adapted <- run_stage(burned, stage = "adapt", iter = 1000, particles = 100, n_cores = 4); sampled <- run_stage(adapted, stage = "sample", iter = 1000, particles = 50, n_cores = 4)
ls()
?usethis::use_data
usethis::use_data(sampled, wgnmks2008Fast, internal=TRUE)
history()
?usethis::use_data_raw
usethis::use_data_raw("
usethis::use_data_raw("wagenmakers2008")
?savehistory
savehistory('wagen.R')
