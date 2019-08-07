# COMPARISON of R code and MATLAB code

# R specific - load libraries and clear workspace
rm(list=ls())
library(rtdists)
library(mvtnorm) ## For the multivariate normal.
library(MASS) ## For matrix inverse.
library(MCMCpack) ## For the inverse Wishart random numbers.
#------------------------------------------------------------------------------------------------
## Load Forstmann et al.'s data.
data=read.csv("data/data.csv",header=FALSE)
names(data)=c("subject","rt","correct","condition")
S=3
#------------------------------------------------------------------------------------------------
## Extract numbers of trials per subject not used in R code
#
#
#
#------------------------------------------------------------------------------------------------
## Grant ability to restart processing
restart=FALSE
cpus=1#5## =1 will give serial operation.
#------------------------------------------------------------------------------------------------
## Set PMwG process "parameters"
n.particles=100
n.iterations=10000
parameter.names=c("b1","b2","b3","A","v1","v2","t0")
n.parameters=length(parameter.names)
# GC: no num_choice
# GC: no burn
# GC: No max iter
#------------------------------------------------------------------------------------------------
## Priors.
prior.mu.mean=rep(0,n.parameters)
prior.mu.sigma2=diag(rep(1,n.parameters))
## Tuning settings for the Gibbs steps.
v.half=2
A.half=1
#------------------------------------------------------------------------------------------------
## Storage for the samples. (GC: And start points)
param.theta.sigma2=array(NA,dim=c(n.parameters,n.parameters,n.iterations),dimnames=list(parameter.names,parameter.names,NULL))
start_points=c(.2,.2,.2,.4,.3,1.3,-2)  # GC: Starting values for parameter mus (population-level parameters)
# GC: We do not do an inverse wishart here, just preallocate NA vals
latent.theta.mu=array(NA,dim=c(n.parameters,S,n.iterations),dimnames=list(parameter.names,NULL,NULL))
# GC: No setting for across-trial variability - seems to be unused in MATLAB code
# GC: MATLAB has duplicate setting for number of random effects - seems to be unused
#------------------------------------------------------------------------------------------------
## GC: R specific steps outside the loop
# GC: Grab a 1 subject slice of the latent array.
param.theta.mu=latent.theta.mu[,1,]
# Make single-iteration-sized versions, for easier reading of code below.
ptm=param.theta.mu[,1]
pts2=param.theta.sigma2[,,1]
# Start points for the population-level parameters only. Hard coded here, just
# for convenience.
ptm[1:n.parameters]=start_points # Weird subscripts maintains the naming.
pts2=diag(rep(.01,n.parameters)) ## Who knows?? # GC: Set some initial values for param.theta.sigma2
pts2.inv=ginv(pts2) ## Because this is calculated near the end of the main loop, needs initialising for iter=1.
## Things I save rather than re-compute inside the loops.
k.half=v.half+n.parameters-1+S
v.shape=(v.half+n.parameters)/2
prior.mu.sigma2.inv=ginv(prior.mu.sigma2) # GC: Inverse of prior
#------------------------------------------------------------------------------------------------
## GC: R version of log-likelihood LBA code is here, this is covered in LBA_MC_v1
# Function to calculate the log-likelihood of data, given random
# effects, or coversely to produce synthetic data from given random
# effects which match shape of "data".
ll=function(x,data,sample=FALSE) {
    x=exp(x)
    if (any(data$rt<x["t0"])) return(-Inf)
    ##b.ind=paste0("b",data$condition)
    ##bs=x["A"]+x[b.ind]
    bs=x["A"]+x[c("b1","b2","b3")][data$condition] ## This is faster than "paste".
    if (sample) {
        out=rLBA(n=nrow(data),A=x["A"],b=bs,t0=x["t0"],mean_v=x[c("v1","v2")],sd_v=c(1,1),distribution="norm",silent=TRUE)
    } else {
        out=dLBA(rt=data$rt,response=data$correct,A=x["A"],b=bs,t0=x["t0"],mean_v=x[c("v1","v2")],sd_v=c(1,1),distribution="norm",silent=TRUE)
        bad=(out<1e-10)|(!is.finite(out))
        out[bad]=1e-10
        out=sum(log(out))
    }
    out
}
#------------------------------------------------------------------------------------------------
## Sample the initial values for the random effects. Algorithm is same
## as for the main resampling down below.
particles=array(dim=c(length(ptm),S))
for (s in 1:S) {
    proposals=rmvnorm(n.particles,ptm,pts2)
    colnames(proposals)=names(ptm) # stripped otherwise.
    lw=apply(proposals,1,ll,data=data[data$subject==s,])
    weight=exp(lw-max(lw))
    particles[,s]=proposals[sample(x=n.particles,size=1,prob=weight),]
}
#------------------------------------------------------------------------------------------------
## GC: More initial values
# GC: No value for tot_param??
## Sample the mixture variables' initial values.
a.half=1/rgamma(n=n.parameters,shape=0.5,scale=1)
# No <count> to keep track of number of subjects meeting criteria to switch
# No <mean_theta> allocated - possibly  related to "single iteration sized bits"
# No <covmat_theta> - similar to above but for covariance matrix
# No <switch_num> - - R code just uses simple naive sampling variant
# No <temp>
# No <i> iterate var
#------------------------------------------------------------------------------------------------
## GC: R code reuseable sampling function
## Paralellisable particle sampling function (cond/orig/real/log).
pm.corl=function(s,data,n.particles,mu,sig2,particles) {
    ## This uses the simplest, and slowest, proposals: mixture of the
    ## the popultion distribution and gaussian around current random effect.
    wmix=0.5
    n1=rbinom(n=1,size=n.particles,prob=wmix)
    if (n1<2) n1=2
    if (n1>(n.particles-2)) n1=n.particles-2 ## These just avoid degenerate arrays.
    proposals1=rmvnorm(n1,mu,sig2)
    proposals2=rmvnorm(n.particles-n1,particles[,s],sig2)
    proposals=rbind(proposals1,proposals2)
    colnames(proposals)=names(mu) # stripped otherwise.
    proposals[1,]=particles[,s] # Put the current particle in slot 1.
    lw=apply(proposals,1,ll,data=data[data$subject==s,]) # Density of data given random effects proposal.
    lp=dmvnorm(x=proposals,mean=mu,sigma=sig2,log=TRUE) # Density of random effects proposal given population-level distribution.
    lm=log(wmix*exp(lp)+(1-wmix)*dmvnorm(x=proposals,mean=particles[,s],sigma=sig2)) # Density of proposals given proposal distribution.
    l=lw+lp-lm # log of importance weights.
    weight=exp(l-max(l))
    proposals[sample(x=n.particles,size=1,prob=weight),]
}
#------------------------------------------------------------------------------------------------
## GC: R code can restart from earlier run
if (restart) {
    cat("\nRestarting from saved run.\n")
    load("data/output/restart.RData")
    ## Check a couple of things.
    if (dim(particles)[1]!=n.parameters) stop("Restart does not match size (params).")
    if (dim(particles)[2]!=S) stop("Restart does not match size (subjects).")
}
#------------------------------------------------------------------------------------------------
## GC: R setup for parallel processing
if (cpus>1) {
    library(snowfall)
    sfInit(parallel=TRUE,cpus=4)
    sfClusterSetupRNG()
    sfLibrary(rtdists)
    ##sfLibrary(mgcv) ## For the multivariate random normal.
    sfLibrary(mvtnorm)
    sfLibrary(MASS) ## For matrix inverse.
    sfLibrary(MCMCpack) ## For the inverse Wishart random numbers.

    sfExportAll(except=c("param.theta.mu","param.theta.sig2","latent.theta.mu"))
}
#------------------------------------------------------------------------------------------------
## GC: Main loop over iterations
for (i in 1:n.iterations) {
    cat("\t",i)
    # Sample population-level parameters. (MATLAB: sample parameter \mu_{\alpha} in Gibbs step, see algorithm 3 (2a) of
    # the paper for the full conditional distribution of the \mu_{\alpha}
    # given the rest of the parameters)
    var_mu=ginv(S*pts2.inv+prior.mu.sigma2.inv)
    mean_mu=as.vector(var_mu%*%(pts2.inv%*%apply(particles,1,sum)))
    chol_var_mu=t(chol(var_mu)) ## t() because I want lower triangle.
    ptm=rmvnorm(1,mean_mu,chol_var_mu%*%t(chol_var_mu))[1,] ## New sample for mu.
    names(ptm)=parameter.names
#------------------------------------------------------------------------------------------------
    ## GC: Sample covariance matrix and sigma2 parameter
    # (MATLAB: sample parameter \Sigma_{\alpha} in Gibbs step, see algorithm 3 (2b)
    # of the paper for the full conditional distribution of the
    # \Sigma_{\alpha} given rest of the parameters)
    # k_half generate earlier outside of loop
    # cov.temp generated using matrix multiplication (MATLAB created more manually)
    theta.temp=particles-ptm
    cov.temp=(theta.temp)%*%(t(theta.temp))
    #
    #
    B.half=2*v.half*diag(1/a.half)+cov.temp
    pts2=riwish(k.half,B.half) ## New sample for sigma.
    pts2.inv=ginv(pts2)
#------------------------------------------------------------------------------------------------
    ## Sample new mixing weights.
    a.half=1/rgamma(n=n.parameters,shape=v.shape,scale=1/(v.half+diag(pts2.inv)+A.half))
    #
    #
    #
    #
    #
#------------------------------------------------------------------------------------------------
## GC: No equivalent for switching to more advanced?? proposal.
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
#------------------------------------------------------------------------------------------------
## Sample new particles for random effects.
    if (cpus>1) {
        tmp=sfLapply(x=1:S,fun=pm.corl,data=data,n.particles=n.particles,mu=ptm,sig2=pts2,particles=particles)
    } else {
        tmp=lapply(X=1:S,FUN=pm.corl,data=data,n.particles=n.particles,mu=ptm,sig2=pts2,particles=particles)
    }
    particles=array(unlist(tmp),dim=dim(particles))
#------------------------------------------------------------------------------------------------
## Store results.
    latent.theta.mu[,,i]=particles
    param.theta.sigma2[,,i]=pts2
    param.theta.mu[,i]=ptm
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
#------------------------------------------------------------------------------------------------
## GC: No counting of unique vals here
    #
    #
    #
    #
    #
#------------------------------------------------------------------------------------------------
## GC: No auto save of data during processing, but end loop
    #
    #
    #
    #
    #
    #
    #
}
#------------------------------------------------------------------------------------------------
## GC: Stop parallel processing and save final dataset
if (cpus>1) sfStop()

save(file="data/output/restart.RData",list=c("pts2.inv","particles","ptm","pts2"))
save.image("data/output/PMwG.RData")
