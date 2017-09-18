function [lambda,mu,delta,FPR,nu,AIC,BIC] = likelihood_preds_m(Data,frames,m,model,FPR,nu)
%{
INPUTS -- 
Data - An N by M binary data matrix, each COLUMN is an independent realisation of the imaging process of an individual fluorophore and is of length M (number of frames). There are N (number of emitters) columns being independenly imaged. The data is binary, so Data(i,j) = 1 implies that the jth fluorophore has been localised in the ith frame and Data(i,j) = 0 implies no localisation. 
frames - Number of frames sampled per second.
m - Number of multiple off states. This program can only estimate parameters for the models with a total number of 1 off state (m=0), 2 off states (m=1) or 3 off states (m=2). If m is not 0,1 or 2, the program will fail.
model - binary vector of length m+2 indicating the model assumed from the data. Since there are m+1 'dark' states 0_0, 0_1, ... 0_m, 1 'on' state 1 and 1 'bleached' state 2, specification of the model will indicate which of the states (in the exact order) 0_0, 0_1, ... 0_m, 1 can access the bleached state 2. For example if m=2 and model = [1 0 1 0], then the bleached state can only be reached by states 0_0 and 0_2, hence the death rates from states 0_1 and 1 are automatically zero and do not need to be estimated. 
FPR - an indication {0,1} if there are false positive observations in the data. If 1, then the false positive rate (FPR) will be estimated and if 0, then it will be automatically set to zero. If this is left blank, then FPR = 1 by default and will be estimated. 
nu - the initial probability mass over the states 0_0, 0_1, ..., 0_m, 1, 2 (at time 0). If KNOWN, this needs to be specified in the form of a vector with the last element zero (since there should be no initial mass on the absorption/bleached state). For example if m=2, nu can be = [0.1 0.2 0 0.7 0]; the sum should be 1. If nu is UNKNOWN input nu=0 and this mass will be estimated. If this is left blank, then nu will be estimated by default. 
%}

%{
OUTPUTS -- 
lambda - a 1 by 2*(m+1) vector of ML estimated transition rates in the order lambda_001, lambda_01, lambda_0102, lambda_011, ..., lambda_0m-10m, lambda_0m-11, lambda_0m1.
mu - a 1 by m+2 vector of ML estimated absorption/bleached rates in the order mu_0, mu_01, .., mu_0m, mu_1. 
delta - the ML estimated noise parameter reliable for false negative observations and setting the lower threshold of photon counts (in time) used to determine whether a localisation has been produced in a frame. 
FPR - the ML estimated false positive rate if not zero. 
nu - either the known or ML estimated initial probability mass of hidden states. 
AIC - reported AIC given the ML estimates of all unknown parameters - can be used for model selection purposes. 
BIC - reported BIC given the ML estimates of all unknown parameters - can be used for model selection purposes. 
%}

if norm(logical(Data)-Data)
	disp('Error: Data matrix must be binary.');
	return 
end 

if m~=0 && m~=1 && m~=2
	disp('Error: Can only estimate parameters for models with multiple off states m=0,1 or 2.');
	return 
end

if length(model)-2 ~= m
	disp('Error: Incorrect specification of model with number of multiple off states m, must be of length m+2.');
	return 
end 
	
if norm(logical(model)-model)
	disp('Error: Incorrect specification of model, model must be a binary vector.');
	return 
end

if nargin == 4
	FPR = 1; %If unspecified, the false positive rate will be a parameter to be estimated.
	nu = 0; %Unspecified, so will be estimated. 
elseif nargin ==5 
	FPR = logical(FPR); %Unless specified as 0, false positive rate will be estimated. 
    nu = 0; %Unspecified, so will be estimated.
elseif nargin == 6
	FPR = logical(FPR); 
	if sum(nu) ~= 0 %Means that nu is known to the experimenter
		if length(nu)-3 ~=m 
			disp('Error: Incorrect specification of initial probability mass, must be of length m+3.');
			return
		end
		
		if sum(nu)~=1 || sum(nu<0)>0
			disp('Error: Incorrect specification of initial probability mass.');
			return
		end 
		
		if nu(end) ~= 0
			disp('Error: Incorrect specification of initial probability mass, cannot start in the bleached state.');
			return
		end 
	end 
elseif nargin < 4
        error('Error: invalid number of input parameters');
end

	
Data = Data'; 
Data = Data(sum(Data')~=0,:); %Filtering out all fluorophores with no observations. 
Delta = 1/frames; %frame rate

if m==0 %Gaining starting parameters. 
	lambda_start = initialisation(Data,Delta,model);
else
	lambdaold = initialisation(Data,Delta,[0 0]);
	lambda_start = ExampleFitDwellTimes_mstate(Data,Delta, m); %used to get starting values for parameters in the multiple off state models.
	lambda_start(end) = lambdaold(2); 
end
 
if sum(nu) ~= 0 %Nu is known to the experimenter. 
	no_nu_params_estimated = 0; 
	if m==0 %Transformation of values to the real line, needed for unconstrained optimisation. 
		zhat = z_trans_m([lambda_start(1:2) lambda_start(3:4) Delta/10 1e-7*FPR],Delta); 
	else 
		zhat = z_trans_m([lambda_start(1:end-1) lambda_start(end).*model Delta/10 1e-7*FPR],Delta); %Transformed initialisation to real line. 
	end 
		other_params = [model 1 FPR];
		myfun = @(z) -eval_loglik_m(Data,z(1:2*(m+1)+length(other_params)),Delta,z_nu_trans_m(nu(1:m+1)),m); %Negative log-likelihood function to be minimised. 
		options = optimset('MaxIter',2e4,'MaxFunEvals',2e4);
        [z,fval,eflag] = fminsearch(myfun,zhat,options); 
		z(isnan(z))=-Inf; %Needed for parameters already known to the experimenter. 
		preds = lam_trans_m(z,Delta); %Transforming back parameter estimates. 
else %Nu is unknown. 
	no_nu_params_estimated = m+1; %m+1 values of nu to be estimated, since nu is of length m+3 and nu(end) = 0.  
	if m==0 
		zhat = [z_trans_m([lambda_start(1:2) lambda_start(3:4) Delta/10 1e-7*FPR],Delta),z_nu_trans_m(0.3)]; %Transformation of values to the real line, needed for unconstrained optimisation. 
	else 
		zhat = [z_trans_m([lambda_start(1:end-1) lambda_start(end).*model Delta/10 1e-7*FPR],Delta),z_nu_trans_m([0.3 1e-3.*ones(1,m)])]; %Transformed initialisation to real line.
	end 
		other_params = [model 1 FPR];
		myfun = @(z) -eval_loglik_m(Data,z(1:2*(m+1)+length(other_params)),Delta,z((2*(m+1)+length(other_params)+1):(3*(m+1)+length(other_params))),m); %Negative log-likelihood function to be minimised. 
		options = optimset('MaxIter',2e4,'MaxFunEvals',2e4);
        [z,fval,eflag] = fminsearch(myfun,zhat,options); 
		z(isnan(z))=-Inf; %Needed for parameters already known to the experimenter. 
		preds = lam_trans_m(z(1:2*(m+1)+length(other_params)),Delta); %Transforming back parameter estimates. 
		nu_preds = nu_trans_m(z(2*(m+1)+length(other_params)+1:end)); %Transforming back nu estimates. 
		nu = [nu_preds 1-sum(nu_preds) 0]; %Estimated nu vector.
end 	

if eflag ~= 1 
    disp('Error: Not converged');
	return
else 
	no_params_estimated = no_nu_params_estimated + length(z(z~=-Inf)); %Total number of estimated parameters. 
    AIC = 2*no_params_estimated + 2*fval; 
    BIC = log(numel(Data))*no_params_estimated + 2*fval; 
end

lambda = preds(1:2*(m+1)); 
mu = preds(2*(m+1)+1:end-2);
delta = preds(end-1); 
FPR = preds(end); 
end 


