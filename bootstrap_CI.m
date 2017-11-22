function [CI_lambda,CI_mu,CI_delta,CI_FPR,CI_nu] = ...
    bootstrap_CI(Data,frames,m,model,FPR,nu,n_boot) 
%USE THIS FUNCTION TO COMPUTE BOOTSTRAP CONFIDENCE INTERVALS  
%{
INPUTS -- 
Data - An M by N binary data matrix, each ROW is an independent realisation of the imaging process of an individual fluorophore and is of length N (number of frames). There are M (number of emitters) columns being independenly imaged. The data is binary, so Data(i,j) = 1 implies that the ith fluorophore has been localised in the jth frame and Data(i,j) = 0 implies no localisation. 
frames - Number of frames sampled per second.
m - Number of multiple off states. This program can only estimate parameters for the models with a total number of 1 off state (m=0), 2 off states (m=1) or 3 off states (m=2). If m is not 0,1 or 2, the program will fail.
model - binary vector of length m+2 indicating the model assumed from the data. Since there are m+1 'dark' states 0_0, 0_1, ... 0_m, 1 'on' state 1 and 1 'bleached' state 2, specification of the model will indicate which of the states (in the exact order) 0_0, 0_1, ... 0_m, 1 can access the bleached state 2. For example if m=2 and model = [1 0 1 0], then the bleached state can only be reached by states 0_0 and 0_2, hence the death rates from states 0_1 and 1 are automatically zero and do not need to be estimated. 
FPR - an indication {0,1} if there are false positive observations in the data. If 1, then the false positive rate (FPR) will be estimated and if 0, then it will be automatically set to zero. If this is left blank, then FPR = 1 by default and will be estimated. 
nu - the initial probability mass over the states 0_0, 0_1, ..., 0_m, 1, 2 (at time 0). If KNOWN, this needs to be specified in the form of a vector with the last element zero (since there should be no initial mass on the absorption/bleached state). For example if m=2, nu can be = [0.1 0.2 0 0.7 0]; the sum should be 1. If nu is UNKNOWN input nu=0 and this mass will be estimated. If this is left blank, then nu will be estimated by default. 
n_boot - number of bootstrap samples to consider. 
%}

%{
OUTPUTS -- 
CI_lambda,CI_mu,CI_delta,CI_FPR,CI_nu will give the 95% bootstrapped CIs
for the paramers lambda, mu, delta, FPR and nu. 
%} 
Delta = 1/frames; 

lambda_boot = zeros(n_boot,2*(m+1));
CI_lambda = zeros(2,2*(m+1)); 

mu_boot = zeros(n_boot,m+2);
CI_mu = zeros(2,m+2); 

delta_boot = zeros(1,n_boot); 
CI_delta = [0 0]; 

FPR_boot = zeros(1,n_boot); 
CI_FPR = [0 0]; 

nu_boot = zeros(n_boot,m+3); 
CI_nu = repmat(nu,2,1); 

for j=1:n_boot 
    X = datasample(Data,length(Data(:,1))); 
    
    if sum(nu)~=0
        if m < 2 
            [lambda_boot(j,:),mu_boot(j,:),delta_boot(j),FPR_boot(j),nu_boot(j,:)] = likelihood_preds_m(X,1/Delta,m,model,FPR,nu);        
        else 
            [lambda_boot(j,:),mu_boot(j,:),delta_boot(j),FPR_boot(j),nu_boot(j,:)] = likelihood_preds_m_manymax(X,1/Delta,m,model,FPR,nu);        
        end 
    else 
        if m < 2 
            [lambda_boot(j,:),mu_boot(j,:),delta_boot(j),FPR_boot(j),nu_boot(j,:)] = likelihood_preds_m(X,1/Delta,m,model,FPR,0);
        else 
            [lambda_boot(j,:),mu_boot(j,:),delta_boot(j),FPR_boot(j),nu_boot(j,:)] = likelihood_preds_m_manymax(X,1/Delta,m,model,FPR,0);
        end 
    end  
end 

Q_lambda = quantile(lambda_boot,linspace(0,1,41)); 
CI_lambda(1,:) = Q_lambda(2,:);  
CI_lambda(2,:) = Q_lambda(end-1,:);

Q_mu = quantile(mu_boot,linspace(0,1,41));
CI_mu(1,:) = Q_mu(2,:);  
CI_mu(2,:) = Q_mu(end-1,:);

Q_delta = quantile(delta_boot,linspace(0,1,41));
CI_delta(1) = Q_delta(2);  
CI_delta(2) = Q_delta(end-1);

Q_FPR = quantile(FPR_boot,linspace(0,1,41)); 
CI_FPR(1) = Q_FPR(2);  
CI_FPR(2) = Q_FPR(end-1);

Q_nu = quantile(nu_boot,linspace(0,1,41)); 
CI_nu(1,:) = Q_nu(2,:);  
CI_nu(2,:) = Q_nu(end-1,:);

end 
