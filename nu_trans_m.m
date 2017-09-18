function nu = nu_trans_m(z_nu)
%This function is used to transform parameters of the initial probability vector nu in the unconstrained likelihood optimisation problem back to their true values. 
nu = z_nu; 
nu(1) = exp(z_nu(1))/((exp(z_nu(1))+1));

if length(z_nu) > 1 
    for i=2:length(z_nu)
        nu(i) = (1-sum(nu(1:i-1)))*exp(z_nu(i))/((exp(z_nu(i))+1)); 
    end 
end 

end 