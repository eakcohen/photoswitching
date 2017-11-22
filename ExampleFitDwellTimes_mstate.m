function [lambda, fval, eflag] = ExampleFitDwellTimes_mstate(X, Delta,m,scal)
%EXAMPLEFITDWELLTIMES Fits dwell time distributions to trace data 
%   Input:
%   data is the 0s and 1s trace data
%   must be an mxn matrix with traces for individual fluorophores in
%   columns such that there are m rows for m fluorophore traces
%   Delta is the exposure time or frame length in seconds
% Rates seem to change with different scal parameter... not sure why 
onTimeMes = [];
offTimeMes = [];
s = size(X);
Xdash = flipud(X.');
ftime = 1:s(1); 

for i=1:s(1) %Calculated when the last fluorophore was observed in each trace. 
     [~,out_last,~] = unique(X(i,:), 'last');
     ftime(i) = out_last(end); 
end 

%ftime(ftime==s(2)) = s(2)-1; 

N = 1:s(1); 
    for j=1:s(1) 
        Y = X(j,1:ftime(j));
        N(j) = length(strfind(Y,[1 0])); 
    end 
N_bar = mean(N); %mean number of 1-0 transitions. 

data = X'; 
for mFluor = 1:s(1) 
    frameStates = data(1:ftime(mFluor),mFluor);
    %list of frames with transitions
    transitions = diff(frameStates);
    StateFrames = diff([0 ;find(transitions)]);
 
    if frameStates(1) % the odd states are on states
        onStates = StateFrames(1:2:end) * Delta;
        offStates = StateFrames(2:2:end) * Delta;
    else % the even states are on states
        offStates = StateFrames(1:2:end) * Delta;
        onStates = StateFrames(2:2:end) * Delta;
    end
     
    onTimeMes = [onTimeMes onStates']; %on times measured. 
    offTimeMes = [offTimeMes offStates']; %off times measured. 
     
end
 
pd_on = fitdist(onTimeMes','Exponential');
pd_off = fitdist(offTimeMes','Exponential');

initial_rate = 1/pd_off.mu; 
if m
	[rates, fval, eflag] = exp_fit_preds(offTimeMes',initial_rate/scal,m); %exponential fitting for multiple off states. 
end 
s1 = 1/pd_on.mu; 
l10 = N_bar*s1/(1+N_bar); 
mu = s1/(1+N_bar); 

%Predicted rates. 
if m==0 
	lambda = [initial_rate l10 mu];
elseif m==1 
    lambda = [rates(2) rates(1) rates(3) l10 mu];
elseif m==2
    lambda = [rates(2) rates(1) rates(4) rates(3) rates(5) l10 mu]; 
end 
end