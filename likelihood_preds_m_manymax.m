function [lambda,mu,delta,FPR,nu,AIC,BIC] = likelihood_preds_m_manymax(Data,frames,m,model,FPR,nu)
%USE THIS FUNCTION IF THE LIKELIHOOD SURFACE IS CLEARLY MULTIMODAL (i.e.
%m=2)
%{
INPUTS -- 
Data - An M by N binary data matrix, each ROW is an independent realisation of the imaging process of an individual fluorophore and is of length N (number of frames). There are M (number of emitters) columns being independenly imaged. The data is binary, so Data(i,j) = 1 implies that the ith fluorophore has been localised in the jth frame and Data(i,j) = 0 implies no localisation. 
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

Data = Data(sum(Data,2)~=0,:); %Filtering out all fluorophores with no observations. 
Delta = 1/frames; %frame rate -- this is incorrect. 
other_params = [model 1 FPR];

if m==0 %Gaining starting parameters. 
lambda_start = initialisation(Data,Delta,model);
else
lambdaold = initialisation(Data,Delta,[logical(sum(model(1:end-1))) logical(sum(model(end)))]);
    diff_scals = [1/100 1/80 1/50 1/40 1/25 1/5 1/3 1/2 1:10]; 
    fval_exp = 1:length(diff_scals); 
    eflag = 1:length(diff_scals); 
    lambda_start = ones(length(diff_scals), 2*m+3); 
    for sc=1:length(diff_scals) 
        [lambda_start(sc,:), fval_exp(sc),eflag(sc)] = ExampleFitDwellTimes_mstate(Data(sum(Data,2)>1,:),Delta, m,diff_scals(sc)); %used to get starting values for parameters in the multiple off state models.
        lambda_start(sc,end-1) = lambdaold(2); 
    end 
    mu_start = [lambdaold(3).*model(1:end-1) lambdaold(4)*model(end)];
    lambda_start = lambda_start(imag(fval_exp)==0 & eflag>0,:);
    fval_exp = fval_exp(imag(fval_exp)==0 & eflag>0); 
end 
%model is of length m + 2

if sum(nu) ~= 0 %Nu is known to the experimenter. 
no_nu_params_estimated = 0; 
if m==0 %Transformation of values to the real line, needed for unconstrained optimisation. 
zhat = z_trans_m([lambda_start(1:2) lambda_start(3:4) Delta/100 1e-7*FPR],Delta); 
else 
        zhat = ones(length(fval_exp),3*m+6); 
        fzs = 1:length(fval_exp); 
        for sc=1:length(fval_exp) 
            zhat(sc,:) = z_trans_m([lambda_start(sc,1:end-1) mu_start Delta/100 1e-7*FPR],Delta); %Transformed initialisation to real line. 
            fzs(sc) = -eval_loglik_m(Data,zhat(sc,:),Delta,z_nu_trans_m(nu(1:m+1)),m);
        end
        zhat = zhat(~isnan(fzs),:);
        fzs = fzs(~isnan(fzs));
        [~, J] = sort(-fzs); 
            if length(J) < 5
                J_inds = J; 
            else 
                J_inds = J(end-4:end); 
            end
         lambda_start(J_inds,:);
         zhat = zhat(J_inds,:); 
end 
myfun = @(z) -eval_loglik_m(Data,z(1:2*(m+1)+length(other_params)),Delta,z_nu_trans_m(nu(1:m+1)),m); %Negative log-likelihood function to be minimised. 
%options = optimset('MaxIter',2e4,'MaxFunEvals',2e4);
problem = createOptimProblem('fminunc','objective',myfun,'x0',zhat(1,:),'options',optimoptions(@fminunc,'Algorithm','quasi-newton','Display','off','MaxIterations',2e4,'MaxFunctionEvaluations',2e4));
ms = MultiStart;
custpts = CustomStartPointSet(zhat);
[z,fval,eflag] = run(ms,problem,custpts); 
z(isnan(z))=-Inf; %Needed for parameters already known to the experimenter. 
preds = lam_trans_m(z,Delta); %Transforming back parameter estimates. 
else %Nu is unknown. 
no_nu_params_estimated = m+1; %m+1 values of nu to be estimated, since nu is of length m+3 and nu(end) = 0.
        nu_start = sum(Data(:,1)==0)/length(Data(:,1));
        if  nu_start == 0 
            nu_start = 1e-6; 
        elseif nu_start == 1 
            nu_start = 1-1e-2;
        end 
if m==0 
    zhat = [z_trans_m([lambda_start(1:2) lambda_start(3:4) Delta/100 1e-7*FPR],Delta),z_nu_trans_m(nu_start)]; %Transformation of values to the real line, needed for unconstrained optimisation. 
else
    zhat = ones(length(fval_exp),3*m+6); s
    fzs = 1:length(fval_exp); 
    for sc=1:length(fval_exp)
        zhat(sc,:) = [z_trans_m([lambda_start(sc,1:end-1) mu_start Delta/100 1e-7*FPR],Delta),z_nu_trans_m([nu_start 1e-6.*ones(1,m)])]; %Transformed initial
        fzs(sc) = -eval_loglik_m(Data,zhat(sc,1:2*(m+1)+length(other_params)),Delta,zhat(sc,(2*(m+1)+length(other_params)+1):(3*(m+1)+length(other_params))),m);
    end
        zhat = zhat(~isnan(fzs),:);
        fzs = fzs(~isnan(fzs));
        [~, J] = sort(-fzs);
            if length(J) < 5
                J_inds = J;
            else
                J_inds = J(end-4:end);
            end
        zhat = zhat(J_inds(end),:);   
end 
myfun = @(z) -eval_loglik_m(Data,z(1:2*(m+1)+length(other_params)),Delta,z((2*(m+1)+length(other_params)+1):(3*(m+1)+length(other_params))),m); %Negative log-likelihood function to be minimised. 
%options = optimoptions(@fminunc,'MaxIterations',2e4,'MaxFunctionEvaluations',2e4);
problem = createOptimProblem('fminunc','objective',myfun,'x0',zhat(1,:),'options', optimoptions(@fminunc,'Algorithm','quasi-newton','Display','off','MaxIterations',2e4,'MaxFunctionEvaluations',2e4));
ms = MultiStart;
custpts = CustomStartPointSet(zhat);
[z,fval,eflag] = run(ms,problem,custpts); 
z(isnan(z))=-Inf; %Needed for parameters already known to the experimenter. 
preds = lam_trans_m(z(1:2*(m+1)+length(other_params)),Delta); %Transforming back parameter estimates. 
nu_preds = nu_trans_m(z(2*(m+1)+length(other_params)+1:end)); %Transforming back nu estimates. 
nu = [nu_preds 1-sum(nu_preds) 0]; %Estimated nu vector.
end

lambda = preds(1:2*(m+1));
mu = preds(2*(m+1)+1:end-2);
delta = preds(end-1);
FPR = preds(end);
AIC = NaN;
BIC = NaN; 

if eflag <= 0
    disp('Error: Not converged');
else 
no_params_estimated = no_nu_params_estimated + length(z(z~=-Inf)); %Total number of estimated parameters. 
    AIC = 2*no_params_estimated + 2*fval; 
    BIC = log(numel(Data))*no_params_estimated + 2*fval; 
end

if bootstrap 
    bootstrap_samples = 
end 

