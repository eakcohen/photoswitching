function [Q_0, Q_1] = emission_delta_endstate_4state(lambda, mu, delta, Delta) 
%This function computes emission matrices when m=1 (2 dark states). It uses the function qij_4state for inverse Laplace transforms and distributions of waiting times. 
%NOTE: We only compute matrices up until n=2 jumps within a time frame, there will be numerical underflow if the number of frames sampled per second is low. 
%m+1 zero states and 1 on state. 
%lambda = (lambda_001 lambda_01 lambda_0102 lambda_011 lambda_0203 lambda_021 ... lambda_0m1 lambda_10
%mu = [mu_0 mu_01 mu_02 ... mu_0m mu_1]
m = length(mu); %m states without the bleached state  

G = zeros(m,m); 
G(1,2) = lambda(1); %l_001
G(1,end) = lambda(2); %l_01
G(end,1) = lambda(end); %l_10

lam_zerom = lambda(3:2:end-1);
lam_ones = [lambda(4:2:end-1) lambda(end-1)];
for i=1:(m-2) 
    G(1+i,end) = lam_ones(i);
    G(1+i,i+2) = lam_zerom(i);
end 

for i=1:m 
    G(i,i) = -(sum(G(i,:))+mu(i));
end 

Q_0 = zeros(m,m+1);
Q_1 = zeros(m,m+1);
Q0_0 = zeros(m,m); 
Q0_1 = Q0_0; 

%%%WHEN X(\DELTA) != 1
for n=0:2
    [Q, endstate] = qij_4state(n,lambda,mu,Delta);
    Q0 = Q(:,[1:m-1 end]); %ends not at one
    Q1 = Q(:,end-1); %ends at 1 
    
    p0 = cdf('Gamma',delta,n,(-G(m,m))^-1)/cdf('Gamma',Delta,n,(-G(m,m))^-1); %starts not at zero 
    p1 = cdf('Gamma',delta,n+1,(-G(m,m))^-1)/cdf('Gamma',Delta,n+1,(-G(m,m))^-1); %starts at one
    
    Q0_0(end,:) = p1*Q0(end,:); 
    Q0_0(1:end-1,:) = p0*Q0(1:end-1,:); 
    Q0_1(end,:) = (1-p1)*Q0(end,:); 
    Q0_1(1:end-1,:) = (1-p0)*Q0(1:end-1,:); 
    
    if n==0 
        Q1_0 = zeros(3,1); 
        Q1_1 = Q1;
    elseif n==1 
        probs_in_state = 1:m-1;
            for i=1:m-1
                probs_in_state(i) = cdf('Gamma',Delta-delta,1,(-G(i,i))^-1)/cdf('Gamma',Delta,1,(-G(i,i))^-1);
            end 
            Q1_1 = [probs_in_state probs_in_state(1)]'.*Q1;
            Q1_0 = (ones(1,m) - [probs_in_state probs_in_state(1)])'.*Q1;
    elseif n==2
			probs = zeros(m,m-1); 
			probs(1,1) = cdf('Gamma',Delta-delta,2,(-G(1,1))^-1)/cdf('Gamma',Delta,2,(-G(1,1))^-1);
		if abs(G(2,2) - G(1,1)) < 1e-4 %To reduce numerical overflow.  
			probs(1,2) = probs(1,1); 
		else %Sum of 2 exponential distributions (closed form). 
			probs(1,2) = 1 - ((-G(2,2)*(1-exp(-G(1,1)*(Delta-delta)))+G(1,1)*(1-exp(-G(2,2)*(Delta-delta))))/(-G(2,2)*(1-exp(-G(1,1)*(Delta)))+G(1,1)*(1-exp(-G(2,2)*(Delta)))));
        end 
			probs(2,2) = probs(1,2); 
			probs(3,:) = probs(1,:);
            
        Q1_1 = sum((probs).*endstate,2);
        Q1_0 = sum((1-probs).*endstate,2);
    end 
    
    Q_0 = Q_0 + [Q0_0(:,1:m-1) Q1_0 Q0_0(:,end)]; %emission matrix when 0 is observed. 
    Q_1 = Q_1 + [Q0_1(:,1:m-1) Q1_1 Q0_1(:,end)]; %emission matrix when 1 is observed. 
end 

Q_0(end+1,:) = zeros(1,m+1);
Q_0(end,end) = 1; 
Q_1(end+1,:) = zeros(1,m+1); 
 
end 


    






