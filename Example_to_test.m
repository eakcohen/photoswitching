%Example to test

load example_data_3state.mat 
load example_data_4state.mat

%3 state simulated under model [0 1] with lambda = [0.3 3], mu = [0 0.1], nu = [0.15 0.85 0], NO false positive observations, with 30 frames sampled per second and a noise threshold of 0.001 seconds (delta).  
%Wish to estimate lambda, mu, nu, delta and the BIC. 

[lambda_3state,mu_3state,delta_3state,~,nu_3state,~,BIC_3state] = likelihood_preds_m(example_data_3state,30,0,[0 1],0); 

%This should give: 
%lambda_3state = [0.2950, 3.0118], mu_3state = [0,0.1013], delta_3state =
%0.0011, nu_3state = [0.1182, 0.8818, 0] and BIC_3state = 5.2286e+04. 

%4 state simulated under model [0 1 0] with lambda = [0.2 0.3 0.1 0.7], mu = [0 0.1 0], nu (KNOWN) = [0 0 1 0], false positive observations with FPR = 1e-6, with 30 frames sampled per second and a noise threshold of 0.001 seconds (delta).  
%Wish to estimate lambda, mu and FPR. 

[lambda_4state,mu_4state,~,FPR_4state] = likelihood_preds_m(example_data_4state,30,1,[0 1 0],1,[0 0 1 0]); 

%This should give: 
%lambda_4state = [0.1654, 0.2804, 0.0834, 0.7125] mu_4state = [0,0.1171,0] and FPR_4state = 1.0541e-06; 

