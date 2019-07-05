function [rnorm_v_kron]=reshape_v(data_response_repmat,rnorm_v1_kron,rnorm_v2_kron)
     
     %this function finds the combination of the mean of the drift rate
     %given the response choice data.
     ind1=data_response_repmat==1;
     ind2=data_response_repmat==2;
      
     rnorm_v_kron(ind1,:)=[rnorm_v1_kron(ind1,1),rnorm_v2_kron(ind1,1)];
     rnorm_v_kron(ind2,:)=[rnorm_v2_kron(ind2,1),rnorm_v1_kron(ind2,1)];

end