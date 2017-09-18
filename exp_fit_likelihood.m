function likelihood = exp_fit_likelihood(lambda,off_times,m)
%This function calculated the log-likelihood for the exponential fitting method when m=1 and 2. 
%lambda = (k_21 k_23 k_31 k_34 k_41) IN THIS FORM 
lambda = exp(lambda); %transformation... 

if m==1 
    l1 = lambda(1) + lambda(2); 
    l2 = lambda(3); %k23 = lambda(2) 
    a1 = 1 + lambda(2)/(l2-l1);
    a2 = -lambda(2)/(l2-l1); 
    
    likelihood = sum(log(a1.*l1.*exp(-l1.*off_times) + a2.*l2.*exp(-l2.*off_times)));
    
elseif m==2 
    l1 = lambda(1) + lambda(2); 
    l2 = lambda(3) + lambda(4); 
    l3 = lambda(5); 

    A = l2-l3; 
    B = l3-l1; 
    C = l1-l2; 
    D = (l2-l3)*l2*l3 + (l3-l1)*l1*l3 + (l1-l2)*l1*l2; 

    a1 = 1 + lambda(2)/(l2-l1) + lambda(2)^2*A/D; 
    a2 = -lambda(2)/(l2-l1) + lambda(2)^2*B/D; 
    a3 = lambda(2)^2*C/D; 

    likelihood = sum(log(a1.*l1.*exp(-l1.*off_times) + a2.*l2.*exp(-l2.*off_times) + a3.*l3.*exp(-l3.*off_times)));
end 
end 

