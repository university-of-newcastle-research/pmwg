%% Load data, no libs, just copied code files----------------------------
%% Log likelihood function in LBA_tpdf, LBA_tcdf files---------------------
%% Parameters unnamed, priors set below-----------------------------------
num_randeffect=7; %the total number of random effects in the LBA model. For Forstmann dataset, we have 7 random effects.
prior_mu_mean=zeros(num_randeffect,1); %the prior for \mu_{\alpha}
prior_mu_sig2=eye(num_randeffect);%the prior for \Sigma_{\alpha}


%% No equivalent of pmwgs option, set hyper parameters here--------------
v_half=2; %the hyperparameters of the prior of \Sigma_{\alpha}
A_half=1; %the hyperparameters of the prior of \Sigma_{\alpha}




%% Initial values for theta_mu/theta_sig2-----------------------------------
param.theta_mu=[0.2695;0.2116;-0.0241;-0.4017;0.2978;1.1230;-1.7128]; %the initial values for parameter \mu
param.theta_sig2=iwishrnd(eye(num_randeffect),20); % the initial values for \Sigma


%% Set up number iters, number burnin etc-----------------------------------
param.sv=1; %set across-trial variability in drift rate to 1 as scaling constant.
param.num_randeffect=7; %total number of random effects (subject-level parameters) in the LBA model. For Forstmann dataset, we have 7 random effects.

num_choice=2; % the number of response options in the task
burn=500; % number of burn in iterations
nit=12000; % we take the last 10000 draws out of 12000 draws
s=100000;% the maximum number of iterations the algorithm will run
%% Initialise the random effects------------------------------------------
[theta_latent]=LBA_MC_v1(data,param,num_subjects,num_trials,num_particles); %obtain initial values of the random effects.

%% Copy of init code-----------------------------------------------------
function [theta_latent]=LBA_MC_v1(data,param,num_subjects,num_trials,num_particles)
% No setup - straight into iterate over subject
parfor j=1:num_subjects
    rnorm=mvnrnd(param.theta_mu',param.theta_sig2,num_particles);
    % Manipulate particles (rnorm) and data into format to run log likelihood function
    lw=real(log(LBA_n1PDF_reparam_real(data_rt_repmat, rnorm_theta_A_kron, rnorm_theta_b_kron, rnorm_theta_v_kron, ones(num_particles*num_trials(j,1),1),rnorm_theta_tau_kron)));
    % Check for real values from LBA function, reshape results
    max_logw=max(real(logw));
    weight=real(exp(logw-max_logw));
    weight=weight./sum(weight);
    % M: Check for neg weights, reshape rnorm_theta etc
    if Nw>0 
        ind=randsample(Nw,1,true,weight);
        theta_latent(j,:)=rnorm_theta(ind,:);
    end
% After looping return theta_latent

%% More setup of vals for gibbs step, storage for conditional step etc --
%% No function for looping, ---------------------------------------------


%% so no checking,-------------------------------------------------------
%% conditional args created later----------------------------------------
while i<=s
%% Gibbs step done inline-------------------------------------------------


%% Last sample/hyper pars already available -----------------------------
    var_mu=inv(num_subjects*inv(param.theta_sig2)+inv(prior_mu_sig2));
    mean_mu=var_mu*(inv(param.theta_sig2)*(sum(theta_latent)'));
    chol_var_mu=chol(var_mu,'lower');
    param.theta_mu=mvnrnd(mean_mu,chol_var_mu*chol_var_mu');
    k_half=v_half+num_randeffect-1+num_subjects;
    cov_temp=zeros(num_randeffect,num_randeffect);
    for j=1:num_subjects
        theta_j=theta_latent(j,:)';
        cov_temp=cov_temp+(theta_j-param.theta_mu')*(theta_j-param.theta_mu')';
    end
    B_half=2*v_half*diag([1./a_half])+cov_temp;
    param.theta_sig2=iwishrnd(B_half,k_half);    
    theta_sig2_inv=inv(param.theta_sig2);

    %sample a_{1},...,a_{7}, see algorithm 3 (2c) of the paper for the full conditional distribution
    for j=1:num_randeffect
        temp_v_half=(v_half+num_randeffect)/2;
        temp_s_half=(v_half*theta_sig2_inv(j,j)+A_half);
        a_half(j,1)=1./random('gam',temp_v_half,1/temp_s_half);
    end
 

%% No equivalent for setting to grab new samples
%% Create efficient (conditional) pars for sampling phase
   % Carlo algorithm
   if sum(count>switch_num)==num_subjects
         for j=1:num_subjects 

             theta=[theta_latent_b1_store(burn:end,j),theta_latent_b2_store(burn:end,j),theta_latent_b3_store(burn:end,j),...
             theta_latent_A_store(burn:end,j),theta_latent_v1_store(burn:end,j),theta_latent_v2_store(burn:end,j),theta_latent_tau_store(burn:end,j),...
             theta_mu_store(burn:end,:),chol_theta_sig2_store1(burn:end,:),chol_theta_sig2_store2(burn:end,:),chol_theta_sig2_store3(burn:end,:),chol_theta_sig2_store4(burn:end,:),chol_theta_sig2_store5(burn:end,:),...
             chol_theta_sig2_store6(burn:end,:),chol_theta_sig2_store7(burn:end,:)];




%% Create covariance matrix and mean for joint random effect and pars----
             covmat_theta(:,:,j)=cov(theta); %computing sample covariance matrix for the joint random effects and parameters \mu_{\alpha} and \Sigma_{\alpha}
             covmat_theta(:,:,j)=topdm(covmat_theta(:,:,j)); %little correction if the covariance matrix for the proposal is not positive definite matrix.
             mean_theta(j,:)=mean(theta); %computing the sample mean for the joint random effects and parameters \mu_{\alpha} and \Sigma_{\alpha}
%% Not done at this point------------------------------------------------

%% Switching num particles etc-------------------------------------------
%% Run conditional monte carlo here--------------------------------------
   [theta_latent]=LBA_CMC_v1(data,param,theta_latent,num_subjects,num_trials,num_particles,mean_theta,covmat_theta,i,count,switch_num);
function [theta_latent,llh]=LBA_CMC_v1(data,param,...
theta_latent,...
num_subjects,num_trials,num_particles,mean_theta,covmat_theta,i,count,switch_num)
%this is the Conditional Monte Carlo algorithm to sample the random effects
%for each subject.
parfor j=1:num_subjects
%% Set scaling parameter, extract last set of random effects (ref par)---        
%% If n.unique good, update mixing pars, num particle to generate and grab conditional mean/var------
% unwind chol_theta_sig2
        xx=[param.theta_mu';chol_theta_sig2_1';chol_theta_sig2_2';chol_theta_sig2_3';...
            chol_theta_sig2_4';chol_theta_sig2_5';chol_theta_sig2_6';chol_theta_sig2_7'];% we need this to compute the mean of the proposal in the sampling stage
        cond_mean=mean_theta(j,1:7)'+covmat_theta(1:7,8:end,j)*((covmat_theta(8:end,8:end,j))\(xx-mean_theta(j,8:end)')); % computing the mean of the proposal in the sampling stage
        cond_var=covmat_theta(1:7,1:7,j)-covmat_theta(1:7,8:end,j)*(covmat_theta(8:end,8:end,j)\covmat_theta(8:end,1:7,j)); % computing the variance of the proposal in the sampling stage
        chol_cond_var=chol(cond_var,'lower');
%condmean = mean_theta + covmat_theta
%
%
%

%% Create proposal particles -------------------------------------------
        rnorm1=cond_mean+chol_cond_var*randn(param.num_randeffect,n1);
        chol_covmat=chol(param.theta_sig2,'lower');
        rnorm2=param.theta_mu'+chol_covmat*randn(param.num_randeffect,n2);
        % Only in 1st two stages (non conditional)
        rnorm1=reference_par'+epsilon.*chol_covmat*randn(7,n1);

        rnorm=[rnorm1,rnorm2];


%% Rearrange arrays ----------------------------------------------------
    rnorm_theta=[rnorm_theta_b1,rnorm_theta_b2,rnorm_theta_b3,rnorm_theta_A,rnorm_theta_v1,rnorm_theta_v2,rnorm_theta_tau];
%% more rearranging to get inn format for :LBA code --------------------
%%computing the log density of the LBA given the particles of random effects---
    lw=real(log(LBA_n1PDF_reparam_real(data_rt_repmat, rnorm_theta_A_kron, rnorm_theta_b_kron, rnorm_theta_v_kron, ones(num_particles*num_trials(j,1),1),rnorm_theta_tau_kron)));
    lw_reshape=reshape(lw,num_trials(j,1),num_particles);
    logw_first=sum(lw_reshape);
    
%% Other densities depending on stage of sampling----------------------------
    if  sum(count>switch_num)==num_subjects
        logw_second=(logmvnpdf(rnorm_theta,param.theta_mu,chol_covmat*chol_covmat'));
        logw_third=log(w_mix.*mvnpdf(rnorm_theta,cond_mean',chol_cond_var*chol_cond_var')+...
            (1-w_mix).*mvnpdf(rnorm_theta,param.theta_mu,chol_covmat*chol_covmat'));
        logw=logw_first'+logw_second'-logw_third;


    else
        logw_second=(logmvnpdf(rnorm_theta,param.theta_mu,chol_covmat*chol_covmat'));
        logw_third=log(w_mix.*mvnpdf(rnorm_theta,reference_par,(epsilon^2).*(chol_covmat*chol_covmat'))+...
            (1-w_mix).*mvnpdf(rnorm_theta,param.theta_mu,chol_covmat*chol_covmat'));
        logw=logw_first'+logw_second'-logw_third;
    end
    
	% mvnpdf(rnorm_theta, ref_par, epsilon^2*(chol_covmat*chol_covmat') 






%%check if there is imaginary number of logw----------------------------
    max_logw=max(real(logw));
    weight=real(exp(logw-max_logw));
    llh_i(j) = max_logw+log(mean(weight)); 
    llh_i(j) = real(llh_i(j)); 	
    weight=weight./sum(weight);
    % More checks for negative weights
%% Sample a winning particle--------------------------------------------
    if Nw>0 
        ind=randsample(Nw,1,true,weight);
        theta_latent(j,:)=rnorm_theta(ind,:);
    end        
    
end
llh=sum(llh_i);
end
   
%% storing the MCMC draws------------------------------------------------
%% Check for number unique values in random effects----------------------
%% Save every N iterations and save final dataset------------------------
