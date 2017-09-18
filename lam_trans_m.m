function lambda = lam_trans_m(z_lambda, Delta)
%This function is used to transform parameters in the unconstrained likelihood optimisation problem back to their true values. 
n = length(z_lambda); %n lambda rates 
lambda = z_lambda; 

for i=1:n-2
    lambda(i) = exp(z_lambda(i)); %lambdas and mus
end 

    lambda(n-1) = Delta*exp(z_lambda(end-1))/((exp(z_lambda(end-1))+1)); %Delta 
    lambda(n) = exp(z_lambda(end))/((exp(z_lambda(end))+1)); %False positive rate 
end 

