
PURPOSE 
The main purpose of this code is to provide Maximum Likelihood (ML) estimates for photo-switching parameters when working with Super-Resolution Imaging Data. 

DATA 
We firstly assume that the data entered is an M by N binary matrix which for simplicity, we will call DATA. 
- Here, each ROW is an independent realisation of the imaging process of an individual fluorophore and is of length N (total number of frames). 
- We further assume that there are M (number of fluorophores) row being INDEPENDENTLY imaged. 
- We need the data to be binary, so DATA(i,j) = 1 implies that the ith fluorophore has been localised in the jth frame and Data(i,j) = 0 implies no localisation of this same fluorophore. 
- Data needs to be in DOUBLE format for the program to work.

THE PROGRAM (PSHMM)

NOTES
- This program will automatically fail if the data entered is not binary. 
- The number of frames sampled per second (which we call 'frames') during the experiment, will need to be provided. 

ITS USES
- Estimate photo-switching parameters (this includes transition rates, absorption/bleached rates and estimates for noise parameters) from a set of candidate models. 
- The basic set of models includes 1 ON/BRIGHT state, m+1 DARK states and 1 permanently-DARK/photo-bleached state. 
- These models are characterised by the term 'm' which when 0 denotes the model of 1 DARK state, when 1 denotes the model of 2 DARK states and when 2 denotes the model of 3 DARK states. 

SPECIFICATIONS
Users can also specify with the binary input 'model' by which of the m+2 transient states (in the order 0_0, 0_1, ... 0_m, 1), the permanently-DARK/photo-bleached state can be reached from. 
- For example, when m=1 and photo-bleaching can only be accessed from the ON/BRIGHT state, the model the user will input would be 'm=1' and 'model = [0 0 1]'. 
- More uses of this can be found in the example_to_test.m file. 

Depending on the type of data inputted, users are also allowed to specify whether there may be false positive observations within the dataset. 
- In this case, users input a '1' for the term 'FPR' and the false positive rate (FPR) can be estimated. 
- If there are no false positive observations, a '0' can be inputted.

Depending on the type of experiment conducted, users are allowed to specify whether the initial probability mass of hidden states (in the order 0_0, 0_1, ..., 0_m, 1, 2 at time 0) is known/unknown to them when gaining the data. 
- In the former case, the program allows for the initial probability mass (in the form of the vector 'nu') to be entered. 
- In the latter, the program allows for the mass 'nu' to be estimated. 
- When entering the mass when it is known, the vector should be of length m+3 and the last element zero as it is assumed that there is no initial probability mass over the permanently-DARK/photo-bleached state and any columns of observations which do not have at least one observation of the fluorophore are automatically removed from the dataset. 
- For example, for experiments conducted when all fluorophores start in the bright state, (e.g. in the case m=2) nu can be entered as nu = [0 0 0 1 0]. 

NAME 
The function that executes the above is called 'likelihood_preds_m.m' and is self-contained. It is not only able to estimate photo-switching parameters for given m, but also is able to report both the AIC and BIC. 
Users are advised (in the case of unknown m) to fit all models (m=0,1,2) and select the model with the lowest AIC/BIC. The BIC is preferred for larger datasets.  

ALTERNATIVE VERSION
For the multimodal version of the program, users are directed to 'likelihood_preds_m_manymax.m' which uses 'fminunc' over a range of starting values and has the above inputs. 

BOOTSTRAP CONFIDENCE INTERVALS
For bootstrapped confidence sets to be reproduced, users are directed to 'bootstrap_CI.m' which has the same inputs as the above functions. 

EXAMPLE
Users are in the first instance, directed to the file 'example_to_test.m' for two examples with two simulated datasets of the usage of this program. 

EXPONENTIAL FITTING 
- This program is in the form 'ExampleFitDwellTimes_mstate.m' and takes in inputs 
- X as the same M by N data as the above
- Delta as 1/frames, where frames are the number of frames sampled per second 
- m, number of multiple off states 
- scal (to address the multimodality) can change the rate predictions. Defaults to 1. 