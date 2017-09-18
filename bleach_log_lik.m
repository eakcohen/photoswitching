function log_lik = bleach_log_lik(ts, lambda,nu,Delta)
%This function evaluates the log-likelihood associated with the absorption rates given times (ts) to absorption. 
lambda = lam_trans_m([lambda -Inf -Inf],Delta); %Transform back from real line. 
lambda(isnan(lambda)) = 0; %For parameters already known to the experimenter. 
T = [-lambda(1)+lambda(3) lambda(1); lambda(2) -(lambda(2)+lambda(4))]; 
t = [lambda(3); lambda(4)]; 

log_lik=0; 
for i=1:length(ts)
    log_lik = log_lik + log(nu * expm(T*ts(i)) * t); %log-likelihood function. 
end 

end 