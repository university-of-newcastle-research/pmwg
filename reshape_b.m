function [rnorm_b_kron]=reshape_b(data_cond_repmat,rnorm_b1_kron,rnorm_b2_kron,rnorm_b3_kron)
    %this function to select appropriate threshold random effects given the
    %conditions of the experiments, for Forstmann we have 3 conditions.
    
    ind1=data_cond_repmat==1;
    ind2=data_cond_repmat==2;
    ind3=data_cond_repmat==3;
    
    rnorm_b_kron(ind1,:)=rnorm_b1_kron(ind1,1);
    rnorm_b_kron(ind2,:)=rnorm_b2_kron(ind2,1);
    rnorm_b_kron(ind3,:)=rnorm_b3_kron(ind3,1);


end