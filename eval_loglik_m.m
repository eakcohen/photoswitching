function likelihood = eval_loglik_m(X,lambda,Delta,nu,m)
%This function evaluates the log-likelihood function for given lambdas, mus and other nusiance paramters. Here lambda is a vector with all transition rates between states, all absorption rates, delta and FPR. 
s = size(X); 
likelihood = 1:s(1);

    lambdas = lam_trans_m(lambda,Delta); %Transform lambda back from real line. 
    lambdas(isnan(lambdas)) = 0; %Needed for parameters already known to the experimenter. 
    nu = nu_trans_m(nu); %Transform nu parameters back from the real line. 
    nu = [nu 1-sum(nu) 0];
    lambda = lambdas(1:2*(m+1)); 
    mu = lambdas(2*(m+1)+1:end-2); %lambdas(end-1) = delta and lambdas(end) = FPR (false positive rate). 

%Need to get emission matrices needed in the log-likelihood from the lambda, mu, nu vectors and from delta (noise parameter). Here delta = lambdas(end-1). 
	if m==0 
		[BP0, BP1] = emission_delta_endstate_3state(lambda, mu, lambdas(end-1),Delta);
	elseif m==1
		[BP0, BP1] = emission_delta_endstate_4state(lambda, mu, lambdas(end-1),Delta);
	elseif m==2
		[BP0, BP1] = emission_delta_endstate_5state(lambda, mu, lambdas(end-1),Delta);
	end
%Update emission matrices to allow for the false positive rate (FPR) which is given by lambdas(end). 
    BP0_FPR = (1-lambdas(end))*BP0; %emission matrix when a 0 is seen. 
    BP1_FPR = BP1 + lambdas(end)*BP0; %emission matrix when a 1 is seen. 

        for i=1:s(1) %For the s(1) independent emitters, apparently this is faster to do in a for loop. 
            likelihood(i) = fwdbkwd_norm_c2(X(i,:),nu,BP0_FPR,BP1_FPR);
        end 
        likelihood = sum(likelihood); %Log-likelihood. 

end 