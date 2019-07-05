## Things that can be useful after the sampling.
## Expect to load some file like the one saved using "save.image"
## at the end of the PMwG.R sampling script.


## Plot chains
par(mfrow=c(4,4),mar=c(2,2,0,0))
for (i in 1:n.parameters) {
    plot(param.theta.mu[i,],type="p",pch=16,cex=.5,col=1)
    matplot(t(latent.theta.mu[i,,]),type="l",pch=16,cex=.5,lty=1)
}
image(latent.theta.mu[1,,-1]==latent.theta.mu[1,,-n.iterations]) ## accept/reject

## Print acceptance rates at the random effect level.
cat("\nAcceptance rates in percent, by subject:\n")
round(100*(1-apply(latent.theta.mu[1,,-1]==latent.theta.mu[1,,-n.iterations],1,mean)))


## Draw some posterior predictive data.
n.posterior=20
burnin=round(n.iterations/2)
synth.data=vector(mode="list",length=n.posterior)
for (i in 1:n.posterior) {
    synth.data[[i]]=data ## To get the right format.
    synth.data$rt=synth.data$correct=rep(NA,nrow(data)) ## Remove real data.
    for (s in 1:S) {
        use=(data$subject==s)
        tmp=ll(x=latent.theta.mu[,s,sample(x=burnin:n.iterations,size=1)],data=data[use,],sample=TRUE)
        synth.data[[i]]$rt[use]=tmp$rt
        synth.data[[i]]$correct[use]=tmp$response
    }
}


## Summary statistics from the data compared to posterior predictive data.
get.acc=function(x) tapply(x$correct==2,x[c("subject","condition")],mean)
get.mrt=function(x) tapply(x$rt,x[c("subject","condition")],mean)
get.qs=function(x) tapply(x$rt,x[c("subject","condition")],quantile,prob=c(.1,.5,.9))
av.acc=get.acc(data)
av.mrt=get.mrt(data)
av.q=array(unlist(get.qs(data)),dim=c(3,dim(av.mrt)),dimnames=c(list(NULL),dimnames(av.mrt)))
pp.acc=pp.mrt=array(NA,dim=c(dim(av.acc),n.posterior))
pp.q=array(NA,dim=c(dim(av.q),n.posterior))
for (i in 1:n.posterior) {
    pp.mrt[,,i]=get.mrt(synth.data[[i]])
    pp.acc[,,i]=get.acc(synth.data[[i]])
    pp.q[,,,i]=unlist(get.qs(synth.data[[i]]))
}

## Plot the quantiles.
tmp.d=apply(av.q,c(1,3),mean)
tmp.pp=apply(pp.q,c(1,3,4),mean)
plot(tmp.d[2,],pch=16,ylim=c(0,1),xlim=c(0,4))
for (i in 1:3) points(y=tmp.pp[i,,],x=rep(1:3,n.posterior),pch=16,col="grey")
arrows(x0=1:3,x1=1:3,y0=tmp.d[1,],y1=tmp.d[3,],angle=90,length=.1,code=3)
points(tmp.d[2,],pch=16)

## Plot the accuracy rates.
tmp.d=apply(av.acc,2,mean)
tmp.pp=apply(pp.acc,c(2,3),mean)
plot(tmp.d,pch=16,ylim=c(0.6,0.9),xlim=c(0,4))
points(y=tmp.pp,x=rep(1:3,n.posterior),pch=16,col="grey")
points(tmp.d,pch=16)

## Correlation matrix from covariance:
tmp=1/sqrt(diag(pts2))
names(tmp)=parameter.names
print(pts2*(tmp%*%t(tmp)),2)

## Compare the prior to posterior, for the marginals only, at popn level.
par(mfrow=c(4,2))
for (i in 1:n.parameters) {
    plot(function(x) dnorm(x,mean=prior.mu.mean[i],sd=prior.mu.sigma2[i,i]^.5),xlim=c(-5,+5))
    lines(density(param.theta.mu[i,]),col="red")
    mtext(side=3,line=-2,at=2,parameter.names[i])
}

## Compare random effects to popn distributions, for marginals only, and
## using only the MAP estimates of each.
tmp1=apply(param.theta.mu[,burnin:n.iterations],1,mean)
tmp2=apply(param.theta.sigma2[,,burnin:n.iterations],1:2,mean)
tmp3=apply(latent.theta.mu[,,burnin:n.iterations],1:2,mean)
par(mfrow=c(4,2))
for (i in 1:n.parameters) {
    plot(function(x) dnorm(x,mean=tmp1[i],sd=tmp2[i,i]^.5),xlim=c(-5,+5))
    points(x=tmp3[i,],y=rep(0,S),pch=4,col="red")
    mtext(side=3,line=-2,at=2,parameter.names[i])
}
par(mfrow=c(4,2))
for (i in 1:n.parameters) {
    plot(x=tmp3[i,],y=data.generating.parameters[i,],pch=16)
    abline(coef=c(0,1))
    mtext(side=3,line=-2,at=2,parameter.names[i])
}

