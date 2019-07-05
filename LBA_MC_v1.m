function [theta_latent]=LBA_MC_v1(data,param,num_subjects,num_trials,num_particles)
% This is the Monte Carlo algorithm for generating initial random effects
% values. This is optional code. We can also set the initial values randomly,
% or use some other methods.

parfor j=1:num_subjects
    % generating the particles for the random effects and adjust the size
    % of the vector of each random effect particles (line 9-24)
    rnorm=mvnrnd(param.theta_mu',param.theta_sig2,num_particles);
    rnorm_theta_b1=rnorm(:,1); 
    rnorm_theta_b2=rnorm(:,2);
    rnorm_theta_b3=rnorm(:,3);
    rnorm_theta_A=rnorm(:,4);
    rnorm_theta_v1=rnorm(:,5);
    rnorm_theta_v2=rnorm(:,6);
    rnorm_theta_tau=rnorm(:,7);
    
    rnorm_theta_b1_kron=kron(rnorm_theta_b1,ones(num_trials(j,1),1));
    rnorm_theta_b2_kron=kron(rnorm_theta_b2,ones(num_trials(j,1),1));
    rnorm_theta_b3_kron=kron(rnorm_theta_b3,ones(num_trials(j,1),1));
    rnorm_theta_A_kron=kron(rnorm_theta_A,ones(num_trials(j,1),1));
    rnorm_theta_v1_kron=kron(rnorm_theta_v1,ones(num_trials(j,1),1));
    rnorm_theta_v2_kron=kron(rnorm_theta_v2,ones(num_trials(j,1),1));
    rnorm_theta_tau_kron=kron(rnorm_theta_tau,ones(num_trials(j,1),1));
    
    %adjust the size of the dataset (line 28-30)
    
    data_response_repmat=repmat(data.response{j,1}(:,1),num_particles,1);
    data_rt_repmat=repmat(data.rt{j,1}(:,1),num_particles,1);
    data_cond_repmat=repmat(data.cond{j,1}(:,1),num_particles,1);  
    
    [rnorm_theta_b_kron]=reshape_b(data_cond_repmat,rnorm_theta_b1_kron,rnorm_theta_b2_kron,rnorm_theta_b3_kron); %choose the threshold particles to match with the conditions of the experiment
    [rnorm_theta_v_kron]=reshape_v(data_response_repmat,rnorm_theta_v1_kron,rnorm_theta_v2_kron); %set the drift rate particles to match with the response
    
    %computing the log density of the LBA given the particles of random effects
    lw=real(log(LBA_n1PDF_reparam_real(data_rt_repmat, rnorm_theta_A_kron, rnorm_theta_b_kron, rnorm_theta_v_kron, ones(num_particles*num_trials(j,1),1),rnorm_theta_tau_kron)));
    lw_reshape=reshape(lw,num_trials(j,1),num_particles);
    logw_first=sum(lw_reshape);
    logw=logw_first';
    
    %check if there is imaginary number of logw
    id=imag(logw)~=0;
    id=1-id;
    id=logical(id);
    logw=logw(id,1); 
    logw=real(logw);

    if sum(isinf(logw))>0 | sum(isnan(logw))>0
     id=isinf(logw) | isnan(logw);
     id=1-id;
     id=logical(id);
     logw=logw(id,1);
    end
    
    max_logw=max(real(logw));
    weight=real(exp(logw-max_logw));
    weight=weight./sum(weight);
    
    if sum(weight<0)>0
        id=weight<0;
        id=1-id;
        id=logical(id);
        weight=weight(id,1);
    end
    Nw=length(weight);
    rnorm_theta=[rnorm_theta_b1,rnorm_theta_b2,rnorm_theta_b3,rnorm_theta_A,rnorm_theta_v1,rnorm_theta_v2,rnorm_theta_tau];
    %choose the random effects according to their weights
    if Nw>0 
        ind=randsample(Nw,1,true,weight);
        theta_latent(j,:)=rnorm_theta(ind,:);
       
    end
        
%----------------------------------------------------------------------------------------------------------------------------------    
    
end

end

