# Demo script

* `create_pmwgs_object(data, par_list, ll_func, priors)`
* [`init_pmwgs_obj(sampler, start points)`](#init_pmwgs)
* [`run_stage(sampler, stage)`](#run_stage)

## init_pmwgs

* get/generate start points for theta\_mu, theta\_sig

```r
for s in subjects:
  particles <- rmvnorm(n.particles, theta_mu, theta_sig)
  # mvnrnd(theta_mu', theta_sig, n), then expand, kronecker tensor product
  lw <- apply(ll_func over particles with subject data)
  weight <- exp(lw - max(lw))
  sample & store in alpha
```

## run_stage

* check and/or set args (mix\_ratio, apply\_fn, [create\_efficient](#create_efficient), parallel etc)
* build sample storage
* general process follows

### R version
```r
for i in iter:
  pars <- new_group_pars(store, x)  # see link below
  tmp <- lapply(new_sample over subjects)
  sm <- unlist(tmp)
  update stage samples w/ pars, sm, ll
  if stage is adapt:
    check_adapted
```

* [`new_group_pars`](#new_group_pars)
* [`new_sample`](#new_sample)
* `update_sampler`

## create_efficient

* Test that necessary conditions met
*

## new_group_pars

* get last sample (`gvi`, `sm`) and hyper parameters
* R version
```r
var_mu <- ginv(n.sub * gvi + prior_sig_inv)  #✔
mean_mu <- var_mu %*% (gvi %*% sum(sm))  #✔
chol_var_mu <- t(chol(var_mu))  #✔
gm <- rmvnorm(1, mean_mu, chol_var_mu %*%t(chol_var_mu))[1,]  #✔
theta_tmp <- sm - gm  #✔
cov_tmp <- theta_tmp %*% (t(theta_tmp))  #✔
B_half <- 2 * v_half * diag(1 / a_half) + cov_tmp  #✔
gv <- riwish(k_half, B_half)  #✔
gvi <- inv(gv)  #✔
a_half <- 1/rgamma(n.pars, v_shape, 1/(v_half + diag(gvi) + A_half))
```
* MATLAB version
```matlab
var_mu=inv(num_subjects*inv(param.theta_sig2)+inv(prior_mu_sig2));
mean_mu=var_mu*(inv(param.theta_sig2)*(sum(theta_latent)'));
chol_var_mu=chol(var_mu,'lower');
param.theta_mu=mvnrnd(mean_mu,chol_var_mu*chol_var_mu');

%sample parameter \Sigma_{\alpha} in Gibbs step, look at the paper for
%the full conditional distribution of the \Sigma_{\alpha} given rest
k_half=v_half+num_randeffect-1+num_subjects;
cov_temp=zeros(num_randeffect,num_randeffect);
for j=1:num_subjects
    theta_j=theta_latent(j,:)';
    cov_temp=cov_temp+(theta_j-param.theta_mu')*(theta_j-param.theta_mu')';
end
B_half=2*v_half*diag([1./a_half])+cov_temp;
param.theta_sig2=iwishrnd(B_half,k_half);    
theta_sig2_inv=inv(param.theta_sig2);

%sample a_{1},...,a_{7}, look at the paper for the full conditional
%distribution
for j=1:num_randeffect
    temp_v_half=(v_half+num_randeffect)/2;
    temp_s_half=(v_half*theta_sig2_inv(j,j)+A_half);
    a_half(j,1)=1./random('gam',temp_v_half,1/temp_s_half);  % random('gam', shape, scale)
end
```
 


## new_sample

* grab efficient pars if necessary, check/extract parameters
* generate new particles with `gen_particles`
* replace first particles with last accepted sample
* Calc log density of data (`lw`) given each random effects proposals with log-likelihood func
  * R version - `apply(proposals, likelihood_func, data for subject)`
  * MATLAB version - `real(log(LBA(data/particles)))` reshaped amd summed
* Calc density of data (`lp`) given population level parameters with multivariate\_normal
  * R version - `dmvnorm(proposals, mu/sig2, log=TRUE)`
  * MATLAB version - `logmvnpdf(particles, theta_mu, chol_covmat*chol_covmat')`
* Calc mixture of densities for each particle type:
  * R version -
    * `prop_density <- dmvnorm(proposals, last_sample_random_effect, sig2)`
    * `eff_density <- dmvnorm(proposals, efficient_mu, efficient_sig2)`
    * `lm <- log( mix[1] * exp(lp) + mix[2] * prop_density + mix[3] * eff_density)`
  * MATLAB version -
    * `log(`
    * `    w1_mix * mvnpdf(data, conditional_mean, chol_cond_var*chol_cond_var') +`
    * `    w2_mix * mvnpdf(particles, theta_mu, chol_covmat*chol_covmat') +`
    * `    w3_mix * mvnpdf(particles, cond_mean_ref', chol_cond_var*chol_cond_var')`
    * `)`
* Combine densities: `lw + lp - lm`

**NOTES:**
* `chol_covmat` is the cholesky factorisation of `theta_sig2` for the lower triangle.
* `chol_cond_var` is the cholesky factorisation of `cond_mean` for the lower triangle.
* `cond_mean_ref` is `theta_latent` for the particular subject, the random\_effects from the previous iteration 
* `cond_mean` is  mean of subject random effects +
