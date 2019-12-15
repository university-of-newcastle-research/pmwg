% estimating the hierarchical LBA model using PMwG method for the Forstmann (2008) dataset
% The LBA specification can be found in the paper: New Estimation
% approaches for the linear Ballistic Accumulator Model
% The data is stored in the matlab file 'LBA_realdata.mat', it has three
% components: 
% data.cond: the conditions of the experiments, we have three conditions in
% the Forstmann data
% data.rt: the response time
% data.response: response = 1 for incorrect and response = 2 for correct.

load('LBA_realdata.mat'); %load the dataset, see an example in the 'LBA_realdata.mat' 
num_subjects=length(data.rt); %number of subjects in the experiments
for j=1:num_subjects
    num_trials(j,1)=length(data.rt{j,1}); %computing the number of trials per subject
end
num_particles=1000; %number of particles in the conditional Monte Carlo algorithm
%parpool(28) %number of processors available to be used.

%initial values of the hyperparameters for lower level parameters
num_randeffect=7; %the total number of random effects in the LBA model. For Forstmann dataset, we have 7 random effects.
prior_mu_mean=zeros(num_randeffect,1); %the prior for \mu_{\alpha}
prior_mu_sig2=eye(num_randeffect);%the prior for \Sigma_{\alpha}
v_half=2; %the hyperparameters of the prior of \Sigma_{\alpha}
A_half=1; %the hyperparameters of the prior of \Sigma_{\alpha}

param.theta_mu=[0.2695;0.2116;-0.0241;-0.4017;0.2978;1.1230;-1.7128]; %the initial values for parameter \mu
param.theta_sig2=iwishrnd(eye(num_randeffect),20); % the initial values for \Sigma
param.sv=1; %set across-trial variability in drift rate to 1 as scaling constant.
param.num_randeffect=7; %total number of random effects (subject-level parameters) in the LBA model. For Forstmann dataset, we have 7 random effects.

num_choice=2; % the number of response options in the task
burn=500; % number of burn in iterations
nit=12000; % we take the last 10000 draws out of 12000 draws
s=100000;% the maximum number of iterations the algorithm will run

[theta_latent]=LBA_MC_v1(data,param,num_subjects,num_trials,num_particles); %obtain initial values of the random effects.

tot_param=42; %total number of parameters for each subjects
a_half=1./random('gam',1/2,1,num_randeffect,1); %initial values for a_{1},...,a_{7}
count=zeros(1,num_subjects);
mean_theta=zeros(num_subjects,tot_param); %allocation for proposal mean for the random effects
covmat_theta=zeros(tot_param,tot_param,num_subjects); %allocation for proposal covariance matrix for the random effects 
switch_num=20;
temp=1;
i=1;

while i<=s
    i
  
    %sample parameter \mu_{\alpha} in Gibbs step, see algorithm 3 (2a) of 
    %the paper for the full conditional distribution of the \mu_{\alpha} 
    %given the rest of the parameters
    var_mu=inv(num_subjects*inv(param.theta_sig2)+inv(prior_mu_sig2));
    mean_mu=var_mu*(inv(param.theta_sig2)*(sum(theta_latent)'));
    chol_var_mu=chol(var_mu,'lower');
    param.theta_mu=mvnrnd(mean_mu,chol_var_mu*chol_var_mu');
    
    %sample parameter \Sigma_{\alpha} in Gibbs step, see algorithm 3 (2b) 
    %of the paper for the full conditional distribution of the 
    %\Sigma_{\alpha} given rest of the parameters
    k_half=v_half+num_randeffect-1+num_subjects;
    cov_temp=zeros(num_randeffect,num_randeffect);
    for j=1:num_subjects
        theta_j=theta_latent(j,:)';
        cov_temp=cov_temp+(theta_j-param.theta_mu')*(theta_j-param.theta_mu')';
    end
    B_half=2*v_half*diag([1./a_half])+cov_temp;
    param.theta_sig2=iwishrnd(B_half,k_half);    
    theta_sig2_inv=inv(param.theta_sig2);

    %sample a_{1},...,a_{7}, see algorithm 3 (2c) of the paper for the 
    %full conditional distribution
    for j=1:num_randeffect
        temp_v_half=(v_half+num_randeffect)/2;
        temp_s_half=(v_half*theta_sig2_inv(j,j)+A_half);
        a_half(j,1)=1./random('gam',temp_v_half,1/temp_s_half);
    end
 
   % training the proposals for conditional Monte
   % Carlo algorithm
   if sum(count>switch_num)==num_subjects
         for j=1:num_subjects 
             theta=[theta_latent_b1_store(burn:end,j),theta_latent_b2_store(burn:end,j),theta_latent_b3_store(burn:end,j),...
             theta_latent_A_store(burn:end,j),theta_latent_v1_store(burn:end,j),theta_latent_v2_store(burn:end,j),theta_latent_tau_store(burn:end,j),...
             theta_mu_store(burn:end,:),chol_theta_sig2_store1(burn:end,:),chol_theta_sig2_store2(burn:end,:),chol_theta_sig2_store3(burn:end,:),chol_theta_sig2_store4(burn:end,:),chol_theta_sig2_store5(burn:end,:),...
             chol_theta_sig2_store6(burn:end,:),chol_theta_sig2_store7(burn:end,:)];
             % in the matrix called theta above, you have to list 
             % (1) all your random effects in the LBA model, in the case of
             % Forstmann, you have \alpha_{b_1}, \alpha_{b_2},
             % \alpha_{b_3}, \alpha_A, \alpha_{v_1}, \alpha_{v_2},
             % \alpha_{tau}, (2) followed by the parameters \mu_{\alpha}, and
             % cholesky factor (lower triangular matrix) of the covariance
             % matrix \Sigma_{\alpha}
             covmat_theta(:,:,j)=cov(theta); %computing sample covariance matrix for the joint random effects and parameters \mu_{\alpha} and \Sigma_{\alpha}
             covmat_theta(:,:,j)=topdm(covmat_theta(:,:,j)); %little correction if the covariance matrix for the proposal is not positive definite matrix.
             mean_theta(j,:)=mean(theta); %computing the sample mean for the joint random effects and parameters \mu_{\alpha} and \Sigma_{\alpha}
         end           
   end
   %if there are at least switch_num unique values of the random effects for each
   %subject, we then switch to better proposal defined in the paper 
   %(equation 12), and reduce the number of particles.
   if sum(count>switch_num)==num_subjects
      num_particles=100; 
      t(temp,1)=i;
      s=t(1,1)+nit;
      temp=temp+1;
   end
   
   %conditional Monte Carlo algorithm to update the random effects, 
   %see algorithm 1 of the paper 
   [theta_latent]=LBA_CMC_v1(data,param,theta_latent,num_subjects,num_trials,num_particles,mean_theta,covmat_theta,i,count,switch_num);
   
    %storing the MCMC draws
    
    % storing the cholesky factor of the covariance matrix \sigma_{\alpha}
    chol_theta_sig2=chol(param.theta_sig2,'lower');
    chol_theta_sig2_store1(i,:)=log(chol_theta_sig2(1,1));
    chol_theta_sig2_store2(i,:)=[chol_theta_sig2(2,1),log(chol_theta_sig2(2,2))];
    chol_theta_sig2_store3(i,:)=[chol_theta_sig2(3,1:2),log(chol_theta_sig2(3,3))];
    chol_theta_sig2_store4(i,:)=[chol_theta_sig2(4,1:3),log(chol_theta_sig2(4,4))];
    chol_theta_sig2_store5(i,:)=[chol_theta_sig2(5,1:4),log(chol_theta_sig2(5,5))];
    chol_theta_sig2_store6(i,:)=[chol_theta_sig2(6,1:5),log(chol_theta_sig2(6,6))];
    chol_theta_sig2_store7(i,:)=[chol_theta_sig2(7,1:6),log(chol_theta_sig2(7,7))];
  
    theta_mu_store(i,:)=param.theta_mu';
    theta_sig2_store1(i,:)=param.theta_sig2(1,:);
    theta_sig2_store2(i,:)=param.theta_sig2(2,2:end);
    theta_sig2_store3(i,:)=param.theta_sig2(3,3:end);
    theta_sig2_store4(i,:)=param.theta_sig2(4,4:end);
    theta_sig2_store5(i,:)=param.theta_sig2(5,5:end);
    theta_sig2_store6(i,:)=param.theta_sig2(6,6:end);
    theta_sig2_store7(i,:)=param.theta_sig2(7,7);
    
    theta_latent_b1_store(i,:)=theta_latent(:,1)';
    theta_latent_b2_store(i,:)=theta_latent(:,2)';
    theta_latent_b3_store(i,:)=theta_latent(:,3)';
    theta_latent_A_store(i,:)=theta_latent(:,4)';
    theta_latent_v1_store(i,:)=theta_latent(:,5)';
    theta_latent_v2_store(i,:)=theta_latent(:,6)';
    theta_latent_tau_store(i,:)=theta_latent(:,7)';
    
    %after burn in, we count the number of unique values in the random
    %effects for each subject. 
    if i>burn
      for j=1:num_subjects
          count(1,j)=size(unique(theta_latent_A_store(burn:end,j)),1);
      end  
    end
    
    %save the output to your directory
     if i==500 | i==1000 | i==2000 | i==3000 | i==4000 | i==5000 | i==6000 | i==7000 | i==8000 | i==9000 | i==10000 | i==20000 | i==30000 | i==40000 | i==50000 | i==60000 | i==70000 | i==80000 | i==90000 | i==100000
        save('/srv/scratch/z3512791/LBA_Forstmann.mat','theta_mu_store','theta_sig2_store1','theta_sig2_store2','theta_sig2_store3','theta_sig2_store4',...
             'theta_sig2_store5','theta_sig2_store6','theta_sig2_store7','theta_latent_b1_store','theta_latent_b2_store','theta_latent_b3_store',...
             'theta_latent_A_store','theta_latent_v1_store','theta_latent_v2_store','theta_latent_tau_store'); 
     end
     i=i+1;   
end
% save the last 10000 draws for further analysis
length_draws=length(theta_latent_A_store);
theta_latent_b1_store=theta_latent_b1_store(length_draws-9999:end,:);
theta_latent_b2_store=theta_latent_b2_store(length_draws-9999:end,:);
theta_latent_b3_store=theta_latent_b3_store(length_draws-9999:end,:);
theta_latent_A_store=theta_latent_A_store(length_draws-9999:end,:);
theta_latent_v1_store=theta_latent_v1_store(length_draws-9999:end,:);
theta_latent_v2_store=theta_latent_v2_store(length_draws-9999:end,:);
theta_latent_tau_store=theta_latent_tau_store(length_draws-9999:end,:);

theta_mu_store=theta_mu_store(length_draws-9999:end,:);
theta_sig2_store1=theta_sig2_store1(length_draws-9999:end,:);
theta_sig2_store2=theta_sig2_store2(length_draws-9999:end,:);
theta_sig2_store3=theta_sig2_store3(length_draws-9999:end,:);
theta_sig2_store4=theta_sig2_store4(length_draws-9999:end,:);
theta_sig2_store5=theta_sig2_store5(length_draws-9999:end,:);
theta_sig2_store6=theta_sig2_store6(length_draws-9999:end,:);
theta_sig2_store7=theta_sig2_store7(length_draws-9999:end,:);

chol_theta_sig2_store1=chol_theta_sig2_store1(length_draws-9999:end,:);
chol_theta_sig2_store2=chol_theta_sig2_store2(length_draws-9999:end,:);
chol_theta_sig2_store3=chol_theta_sig2_store3(length_draws-9999:end,:);
chol_theta_sig2_store4=chol_theta_sig2_store4(length_draws-9999:end,:);
chol_theta_sig2_store5=chol_theta_sig2_store5(length_draws-9999:end,:);
chol_theta_sig2_store6=chol_theta_sig2_store6(length_draws-9999:end,:);
chol_theta_sig2_store7=chol_theta_sig2_store7(length_draws-9999:end,:);

%save the output to your directory. Update directory path as required.
save('/YourDirectoryHere/LBA_Forstmann.mat','theta_mu_store','theta_sig2_store1','theta_sig2_store2','theta_sig2_store3','theta_sig2_store4',...
             'theta_sig2_store5','theta_sig2_store6','theta_sig2_store7','theta_latent_b1_store','theta_latent_b2_store','theta_latent_b3_store',...
             'theta_latent_A_store','theta_latent_v1_store','theta_latent_v2_store','theta_latent_tau_store',...
             'chol_theta_sig2_store1','chol_theta_sig2_store2','chol_theta_sig2_store3','chol_theta_sig2_store4',...
             'chol_theta_sig2_store5','chol_theta_sig2_store6','chol_theta_sig2_store7'); 
     
