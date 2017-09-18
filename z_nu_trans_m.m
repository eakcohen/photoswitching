function z_nu = z_nu_trans_m(nu)
%This function is used to transform parameters of the initial probability vector nu in the unconstrained likelihood optimisation to the real line. 
z_nu = nu; 
z_nu(1) = log(nu(1)/(1-nu(1)));

if length(nu) > 1 
    for i=2:length(nu)
        z_nu(i) = log(nu(i)/(1-sum(nu(1:i)))); %Need sum(nu) = 1. 
    end 
end    

end 