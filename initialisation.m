function lambdahat = initialisation(X,Delta,model) 
%This can be used to find initial estimates in the m=0 3 state case. 
s = size(X); 
ftime = 1:s(1); 

for i=1:s(1) %Find the last time a fluorophore has been seen in each trace. 
     [~,out_last,~] = unique(X(i,:), 'last');
     ftime(i) = out_last(end); 
end 

n00 = 1:s(1); 
n01 = 1:s(1); 
n10 = 1:s(1); 
n11 = 1:s(1); 

%Counting the number of transitions
for i=1:s(1)
        n00(i) = sum(diff(find(X(i,1:ftime(i))==0))==1);
        n01(i) = sum(diff(X(i,1:ftime(i)))>0);
        n10(i) = sum(diff(X(i,1:ftime(i)))<0);
        n11(i) = sum(diff(find(X(i,1:ftime(i))==1))==1);
end 

%Crude MLE estimators for starting values. 
phat = sum(n00(n00~=0))/(sum(n00(n00~=0))+sum(n01(n01~=0))); 
qhat = sum(n10(n10~=0))/(sum(n10(n10~=0))+sum(n11(n11~=0))); 

if phat == 0 
    xbar = mean((s(2)-ftime)*Delta);
    lambdahat1 = 1/(xbar*(1+(1/mn))); 
    lambdahat0 = ((lambdahat1*Delta/qhat) - qhat +1)/Delta; 
else 
lambdahat0 = -log(phat)/Delta; 
lambdahat1 = qhat*lambdahat0/((1-exp(-lambdahat0*Delta))*(exp(-lambdahat0*Delta)-qhat));
end

muhat0 = 0; 
muhat1 = 0; 

h= 0;
for i=1:s(1) 
    h = h + sum(X(i,1:s(2)-ftime(i)))/s(1);
end 

mu_start0 = 1/(lambdahat0*(h*Delta+1)); %starting values of absorption rates. 
mu_start1 = 1/(lambdahat1*(h*Delta+1));
zhat = z_trans_m([lambdahat0 lambdahat1 [mu_start0 mu_start1].*model 0 0],Delta); %Tranform parameters to real line. 

ts=(s(2)-ftime).*Delta; %approximate times to bleached.
if sum(model) > 0 
	myfun = @(z) -bleach_log_lik(ts, [zhat(1) zhat(2) z(1) z(2)],[0.2 0.8],Delta); %Function to minimise. 
	z = fminsearch(myfun,[zhat(3) zhat(4)]);
	z(isnan(z)) = -Inf; 
	muhat0 = exp(z(1)); %Transform parameters back. 
    muhat1 = exp(z(2)); 
end 

 lambdahat = [lambdahat0 lambdahat1 muhat0 muhat1]; 
end 