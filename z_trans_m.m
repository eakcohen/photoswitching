function z = z_trans_m(lambda,Delta)
%This function is used to transform parameters in the unconstrained likelihood optimisation to the real line. 
n = length(lambda); %n lambda rates 
z = lambda; 

for i=1:n-2
    z(i) = log(lambda(i)); %lambdas and mus
end 
    z(n-1) = -log((Delta/(lambda(end-1)))-1); %delta
    z(n) = -log((1/(lambda(end)))-1); %false positive rate
    
end 
