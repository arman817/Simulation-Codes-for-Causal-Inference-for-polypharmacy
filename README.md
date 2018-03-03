## Simulation-Codes-for-Causal-Inference-for-polypharmacy

This repository contains the codes for the paper on "Causal Inference for polypharmacy: Propensity score estimation with multiple concurrent medications" by ARMAN ALAM SIDDIQUE, MIREILLE E. SCHNITZER, ASMA BAHAMYIROU, GUANBO WANG and ANDREA BENEDETTI

Code1.R contains the code for the parameter estimation while dealing with 4 medications and Code2.R and Code3.R contains codes for the parameter estimation for dealing with 8 medications (Case when some of these combinations of the medications aren't observed in the study)

In order to change the sample size from 500 to 1000 increase the value of n from 500 to 1000 in Code1.R

Code2.R corresponds to the estimation of the target parameter using Softmax Regression and Support Vector Machines for the entire data sample

Code3.R refers to the estimation of then target parameter using the Generalized Boosted Method on the entire data sample and using Softmax Regression and Support Vector Machines on the subsetted data.
