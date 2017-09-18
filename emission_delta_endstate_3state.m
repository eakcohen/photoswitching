function [Q_0, Q_1] = emission_delta_endstate_3state(lambda, mu, delta, Delta) 
%This function computes emission probabilities when m=0. It uses inverse laplace transforms from the function qij_3state and distribution of waiting times. 

lambda = [lambda mu]; 
B = zeros(3,2,3); 
if abs(sum(lambda(1)+lambda(3))-sum(lambda(2)+lambda(4)))>1e-3 %When sigma_0 is close to sigma_1. This is used to reduce numerical overflow. 
    q000 = exp(-(lambda(1)+lambda(3))*Delta);
    q020 = lambda(3)*(1-exp(-(lambda(1)+lambda(3))*Delta))/(lambda(1)+lambda(3));
    q100 = (lambda(2)*exp(-(lambda(1)+lambda(3))*Delta)/(lambda(2)+lambda(4)-lambda(1)-lambda(3))) + lambda(2)*exp(-(lambda(2)+lambda(4))*Delta)/(-lambda(2)-lambda(4)+lambda(1)+lambda(3));
    q110 = exp(-(lambda(2)+lambda(4))*Delta); 
    q120 = (lambda(2)*lambda(3)*(1-exp(-(lambda(1)+lambda(3))*Delta))/((lambda(2)+lambda(4)-lambda(1)-lambda(3))*(lambda(1)+lambda(3)))) + (lambda(2)*lambda(3)*(1-exp(-(lambda(2)+lambda(4))*Delta))/((-lambda(2)-lambda(4)+lambda(1)+lambda(3))*(lambda(2)+lambda(4))))+lambda(4)*(1-exp(-(lambda(2)+lambda(4))*Delta))/(lambda(2)+lambda(4));
else 
    q000 = qij_3state(0,0,0,lambda,Delta);
    q020 = qij_3state(0,2,0,lambda,Delta);
    q100 = qij_3state(1,0,0,lambda,Delta);
    q110 = qij_3state(1,1,0,lambda,Delta);
    q120 = qij_3state(1,2,0,lambda,Delta);
end 

%STARTS AT ZERO 
    b000 = q000;
    b001 = 0; 
    diff = 1; 
    i=1; 
        while diff>1e-4 
             p1 = cdf('Gamma',delta,i,(lambda(2)+lambda(4))^-1)/cdf('Gamma',Delta,i,(lambda(2)+lambda(4))^-1);
             add0 = p1*qij_3state(0,0,i,lambda,Delta);
             add1 = (1-p1)*qij_3state(0,0,i,lambda,Delta);
             diff = abs(add0); 
             b000 = b000 + add0; 
             b001 = b001 + add1; 
             i=i+1; 
        end 
        B(1,1,1) = b000; 
        B(1,2,1) = b001; 
        
    b010 = 0; 
    b011 = 0; 
    diff = 1; 
    i=1; 
        while diff>1e-4 
             p2 = 1-(cdf('Gamma',Delta-delta,i,(lambda(1)+lambda(3))^-1))/cdf('Gamma',Delta,i,(lambda(1)+lambda(3))^-1);
             add0 = p2*qij_3state(0,1,i,lambda,Delta);
             add1 = (1-p2)*qij_3state(0,1,i,lambda,Delta); 
             diff = abs(add0); 
             b010 = b010 + add0;
             b011 = b011 + add1; 
             i=i+1; 
        end 
        B(1,1,2) = b010; 
        B(1,2,2) = b011; 
        
    b020 = q020; 
    b021 = 0; 
    diff = 1; 
    i=1; 
        while diff>1e-4 
             p1 = cdf('Gamma',delta,i,(lambda(2)+lambda(4))^-1)/cdf('Gamma',Delta,i,(lambda(2)+lambda(4))^-1);
             add0 = p1*qij_3state(0,2,i,lambda,Delta);
             add1 = (1-p1)*qij_3state(0,2,i,lambda,Delta); 
             diff = abs(add0); 
             b020 = b020 + add0;
             b021 = b021 + add1; 
             i=i+1; 
        end 
        B(1,1,3) = b020; 
        B(1,2,3) = b021; 
     
    %STARTS AT ONE
    b100 = (1-exp(-(lambda(2)+lambda(4))*delta))*(q100)/(1-exp(-(lambda(2)+lambda(4))*Delta)); 
    b101 = q100*(1 - ((1-exp(-(lambda(2)+lambda(4))*delta))/(1-exp(-(lambda(2)+lambda(4))*Delta)))); 
    diff = 1; 
    i=1; 
        while diff>1e-4 
             p1 = cdf('Gamma',delta,i+1,(lambda(2)+lambda(4))^-1)/cdf('Gamma',Delta,i+1,(lambda(2)+lambda(4))^-1);
             add0 = p1*qij_3state(1,0,i,lambda,Delta);
             add1 = (1-p1)*qij_3state(1,0,i,lambda,Delta);
             diff = abs(add0); 
             b100 = b100 + add0; 
             b101 = b101 + add1; 
             i=i+1; 
        end 
        B(2,1,1) = b100; 
        B(2,2,1) = b101; 
    
    b110 = 0; 
    b111 = q110; 
    diff = 1; 
    i=1; 
        while diff>1e-4 
             p2 = (cdf('Gamma',Delta,i,(lambda(1)+lambda(3))^-1)-cdf('Gamma',Delta-delta,i,(lambda(1)+lambda(3))^-1))/cdf('Gamma',Delta,i,(lambda(1)+lambda(3))^-1);
             add0 = p2*qij_3state(1,1,i,lambda,Delta);
             add1 = (1-p2)*qij_3state(1,1,i,lambda,Delta);
             diff = abs(add0); 
             b110 = b110 + add0;
             b111 = b111 + add1; 
             i=i+1; 
        end 
        B(2,1,2) = b110; 
        B(2,2,2) = b111; 
        
    b120 = (1-exp(-(lambda(2)+lambda(4))*delta))*(q120)/(1-exp(-(lambda(2)+lambda(4))*Delta)); 
    b121 = q120*(1-((1-exp(-(lambda(2)+lambda(4))*delta))/(1-exp(-(lambda(2)+lambda(4))*Delta)))); 
    diff = 1; 
    i=1; 
        while diff>1e-4 
             p1 = cdf('Gamma',delta,i+1,(lambda(2)+lambda(4))^-1)/cdf('Gamma',Delta,i+1,(lambda(2)+lambda(4))^-1);
             add0 = p1*qij_3state(1,2,i,lambda,Delta);
             add1 = (1-p1)*qij_3state(1,2,i,lambda,Delta);
             diff = abs(add0); 
             b120 = b120 + add0;
             b121 = b121 + add1; 
             i=i+1; 
        end 
        B(2,1,3) = b120; 
        B(2,2,3) = b121; 
        
   %STARTS AT DEAD 
   B(3,1,3) = 1; 
   
    Q_0 = [B(:,1,1) B(:,1,2) B(:,1,3)]; 
    Q_1 = [B(:,2,1) B(:,2,2) B(:,2,3)];

   if lambda(3)==0 && lambda(4)==0
       Q_0 = [B(1:2,1,1) B(1:2,1,2)];
	   Q_1 = [B(1:2,2,1) B(1:2,2,2)];
   end 
end 
        
    
    