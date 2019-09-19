rm(list=ls())
setwd("~/Documents/Research/Modelling Project/Scripts")
library(rtdists)
##library(mgcv) ## For the multivariate random normal.
library(mvtnorm) ## For the multivariate normal.
library(MASS) ## For matrix inverse.
library(MCMCpack) ## For the inverse Wishart random numbers.

restart=FALSE
cpus=4 ## =1 will give serial operation.
n.particles=200
n.iterations=500
burnin = 250  # REILLY
proposal_parameters = NULL  # REILLY

## Load Forstmann et al.'s data.
data=read.csv("data.csv",header=FALSE)
names(data)=c("subject","rt","correct","condition")
S=19

parameter.names=c("b1","b2","b3","A","v1","v2","t0")
n.parameters=length(parameter.names)

## Storage for the samples.
latent.theta.mu=array(NA,dim=c(n.parameters,S,n.iterations),dimnames=list(parameter.names,NULL,NULL))
param.theta.mu=latent.theta.mu[,1,]
param.theta.sigma2=array(NA,dim=c(n.parameters,n.parameters,n.iterations),dimnames=list(parameter.names,parameter.names,NULL))
proposal_means=array(dim = c(n.parameters,S))  # REILLY
proposal_sigmas=array(dim = c(n.parameters,n.parameters,S))  # REILLY
proposal_parameters=array(NA,dim=c(S,n.parameters,n.parameters,n.parameters),dimnames=list(NULL,parameter.names,parameter.names,parameter.names))  # REILLY
## Make single-iteration-sized versions, for easier reading of code below.
ptm=param.theta.mu[,1]
pts2=param.theta.sigma2[,,1]

## Start points for the population-level parameters only. Hard coded here, just
## for convenience.
ptm[1:n.parameters]=c(.2,.2,.2,.4,.3,1.3,-2) # Weird subscripts maintains the naming.
pts2=diag(rep(.01,n.parameters)) ## Who knows??
pts2.inv=ginv(pts2) ## Because this is calculated near the end of the main loop, needs initialising for iter=1.

## Priors.
prior.mu.mean=rep(0,n.parameters)
prior.mu.sigma2=diag(rep(1,n.parameters))

## Tuning settings for the Gibbs steps.
v.half=2
A.half=1

## Things I save rather than re-compute inside the loops.
k.half=v.half+n.parameters-1+S
v.shape=(v.half+n.parameters)/2
prior.mu.sigma2.inv=ginv(prior.mu.sigma2)

## Function to calculate the log-likelihood of data, given random
## effects, or coversely to produce synthetic data from given random
## effects which match shape of "data".
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

## Sample the mixture variables' initial values.
a.half=1/rgamma(n=n.parameters,shape=0.5,scale=1)


# get.proposals=function(){
# all_samples = rbind(alpha_samples,mu_samples,sigma_samples_unwound)
# mu_tilde = apply(all_samples,1,mean)
# sigma_tilde = var(t(all_samples))
# mu_proposal = rmvnorm(1,mu,sigma)
# sigma_proposal = sigma_samples[,,1] 
# proposal_parameters = condMVN(mean=mu_tilde,sigma=sigma_tilde,dependent.ind=1:length(mu),given.ind=(length(mu)+1):length(mu_tilde),X.given=c(mu_proposal,unwind(sigma_proposal)))
# 
# }


## Paralellisable particle sampling function (cond/orig/real/log).
pm.corl=function(s,data,n.particles,mu,sig2,particles) {
  if (is.na(proposal_means[1,])==TRUE) {  # REILLY
    ## This uses the simplest, and slowest, proposals: mixture of the
    ## the popultion distribution and gaussian around current random effect.
    wmix=0.5
    n1=rbinom(n=1,size=n.particles,prob=wmix)
    if (n1<2) n1=2
    if (n1>(n.particles-2)) n1=n.particles-2 ## These just avoid degenerate arrays.
    proposals1=rmvnorm(n1,mu,sig2)
    proposals2=rmvnorm(n.particles-n1,particles[,s],sig2/2)  # REILLY
    proposals=rbind(proposals1,proposals2)
  }  # REILLY START
  else{
    wmix=c(0.01,0.29,0.7)
    n1=rbinom(n=2,size=n.particles,prob=wmix)
    n2=n1[2]
    n1=n1[1]
    if (n1<2) n1=2
    if (n2<2) n2=2
    if (n1+n2>(n.particles-2)) n1=n.particles-2 ## These just avoid degenerate arrays.
    proposals1=rmvnorm(n1,mu,sig2)
    proposals2=rmvnorm(n2,particles[,s],sig2/2)
    proposals3=rmvnorm(n.particles-(n1+n2),mean=proposal_means[,s],sigma=proposal_sigmas[,,s])
    proposals=rbind(proposals1,proposals2,proposals3)
  }  # REILLY END
    colnames(proposals)=names(mu) # stripped otherwise.
    proposals[1,]=particles[,s] # Put the current particle in slot 1.
    lw=apply(proposals,1,ll,data=data[data$subject==s,]) # Density of data given random effects proposal.
    lp=dmvnorm(x=proposals,mean=mu,sigma=sig2,log=TRUE) # Density of random effects proposal given population-level distribution.
    lm=log(wmix*exp(lp)+(1-wmix)*dmvnorm(x=proposals,mean=particles[,s],sigma=sig2/2)) # Density of proposals given proposal distribution.  # REILLY
    l=lw+lp-lm # log of importance weights.
    weight=exp(l-max(l))
    proposals[sample(x=n.particles,size=1,prob=weight),]
  }  # REILLY

if (restart) {
    cat("\nRestarting from saved run.\n")
    load("restart.RData")
    ## Check a couple of things.
    if (dim(particles)[1]!=n.parameters) stop("Restart does not match size (params).")
    if (dim(particles)[2]!=S) stop("Restart does not match size (subjects).")
}

  # REILLY START
unwind=function(x,reverse=FALSE) {
  ## Takes a variance matrix and unwinds to vector
  ## via Cholesky then log. Or puts it back if
  ## reverse==TRUE.
  if (reverse) {
    ##        if ((n*n+n)!=2*length(x)) stop("Wrong sizes in unwind.")
    n=sqrt(2*length(x)+0.25)-0.5 ## Dim of matrix.
    out=array(0,dim=c(n,n))
    out[lower.tri(out,diag=TRUE)]=x
    diag(out)=exp(diag(out))
    out=out%*%t(out)
  } else {
    y=t(chol(x))
    diag(y)=log(diag(y))
    out=y[lower.tri(y,diag=TRUE)]
  }
  return(out)
}

source("condMNV.R")
  # REILLY END
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

for (i in 1:n.iterations) {
  if (i > burnin){  # REILLY START
    if (i%%20==1) {
      for (s in 1:S){
        pts2.unwound = apply(param.theta.sigma2[,,1:i-1],3,unwind)
        all_samples = rbind(latent.theta.mu[,s,1:i-1],param.theta.mu[,1:i-1],pts2.unwound)
        mu_tilde = apply(all_samples,1,mean)
        sigma_tilde = var(t(all_samples))
        condmvn = condMVN(mean=mu_tilde,sigma=sigma_tilde,dependent.ind=1:length(ptm),given.ind=(length(ptm)+1):length(mu_tilde),X.given=c(ptm,unwind(pts2)))
        proposal_means[,s] = condmvn$condMean
        proposal_sigmas[,,s] = condmvn$condVar
      }
    }
  }  # REILLY END
    cat("\t",i)
    ## Sample population-level parameters.
    var_mu=ginv(S*pts2.inv+prior.mu.sigma2.inv)
    mean_mu=as.vector(var_mu%*%(pts2.inv%*%apply(particles,1,sum)))
    chol_var_mu=t(chol(var_mu)) ## t() because I want lower triangle.
    ptm=rmvnorm(1,mean_mu,chol_var_mu%*%t(chol_var_mu))[1,] ## New sample for mu.
    names(ptm)=parameter.names

    theta.temp=particles-ptm
    ##cov.temp=array(0,dim=c(n.parameters,n.parameters))
    ##for (j in 1:S) cov.temp=cov.temp+(theta.temp[,j])%*%t(theta.temp[,j])
    cov.temp=(theta.temp)%*%(t(theta.temp))
    B.half=2*v.half*diag(1/a.half)+cov.temp
    pts2=riwish(k.half,B.half) ## New sample for sigma.
    pts2.inv=ginv(pts2)

    ## Sample new mixing weights.
    a.half=1/rgamma(n=n.parameters,shape=v.shape,scale=1/(v.half+diag(pts2.inv)+A.half))

    ## Sample new particles for random effects.
    if (cpus>1) {
        tmp=sfLapply(x=1:S,fun=pm.corl,data=data,n.particles=n.particles,mu=ptm,sig2=pts2,particles=particles)
    } else {
        tmp=lapply(X=1:S,FUN=pm.corl,data=data,n.particles=n.particles,mu=ptm,sig2=pts2,particles=particles)
    }
    particles=array(unlist(tmp),dim=dim(particles))

    ## Store results.
    latent.theta.mu[,,i]=particles
    param.theta.sigma2[,,i]=pts2
    param.theta.mu[,i]=ptm
}
if (cpus>1) sfStop()

# for (s in 1:S){
# pts2.unwound = apply(param.theta.sigma2,3,unwind)
# all_samples = rbind(latent.theta.mu[,s,],param.theta.mu,pts2.unwound)
# mu_tilde = apply(all_samples,1,mean)
# sigma_tilde = var(t(all_samples))
# condmvn = condMVN(mean=mu_tilde,sigma=sigma_tilde,dependent.ind=1:length(ptm),given.ind=(length(ptm)+1):length(mu_tilde),X.given=c(ptm,unwind(pts2)))
# proposal_means[,s] = condmvn$condMean
# proposal_sigmas[,,s] = condmvn$condVar
# }

save(file="restart.RData",list=c("pts2.inv","particles","ptm","pts2"))
save.image("PMwG.RData")


#####Checking particle movement#####

#  REILLY START
apply(latent.theta.mu[1,,-1]==latent.theta.mu[1,,-250],1,sum)
apply(latent.theta.mu[1,,-251]==latent.theta.mu[1,,-500],1,sum)

matplot(t(latent.theta.mu[,15,]),type="l")

dev.off()
par(mfcol=c(5,4),mar=c(2,2,1,1))
for (i in 1:s) {
  matplot(t(latent.theta.mu[,i,]),type="l") 
}
# REILLY END
