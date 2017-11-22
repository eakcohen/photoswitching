function likelihood = exp_fit_likelihood(rates,off_times,m)
%This function calculated the log-likelihood for the exponential fitting method when m=1 and 2. 
%rates = (k_21 k_23 k_31 k_34 k_41) IN THIS FORM 
rates = exp(rates); %transformation... 

if m==1 
    l1 = rates(1) + rates(2); 
    l2 = rates(3); %k23 = rates(2) 
    a1 = 1 + rates(2)/(l2-l1);
    a2 = -rates(2)/(l2-l1); 
    
    likelihood = sum(log(a1.*l1.*exp(-l1.*off_times) + a2.*l2.*exp(-l2.*off_times)));
    
elseif m==2 
    l1 = rates(1) + rates(2); 
    l2 = rates(3) + rates(4); 
    l3 = rates(5); 

    A = l2-l3; 
    B = l3-l1; 
    C = l1-l2; 
    D = (l2-l3)*l2*l3 + (l3-l1)*l1*l3 + (l1-l2)*l1*l2; 

    a1 = 1 + rates(2)/(l2-l1) + rates(2)*rates(4)*A/D; 
    a2 = rates(2)/(l1-l2) + rates(2)*rates(4)*B/D; 
    a3 = rates(2)*rates(4)*C/D; 

    likelihood = sum(log(a1.*l1.*exp(-l1.*off_times) + a2.*l2.*exp(-l2.*off_times) + a3.*l3.*exp(-l3.*off_times)));
%     if likelihood < 0 
%         likelihood = 0; 
%     end 
end 
end 

