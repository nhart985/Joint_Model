Tables1_2_Simulations.R: Example simulation code used for Tables 1 and 2 of the paper. The sample size "n" was changed to create each setting. Iterations were run in parallel on a cluster. 

Table3_Simulations.R: Example simulation code used for Table 3 of the paper. The sample size "n" was changed to create each setting. Iterations were run in parallel on a cluster. 

Table4_Simulations.R: Example simulation code used for Table 4 of the paper. Iterations were run in parallel on a cluster. 

Simulations.slurm: Example BASH script to run simulations on cluster.

Combine.slurm: Example BASH script to combine simulation results on cluster.

Combine.R: Code to combine simulation results across iterations (from the parallel comptuing output)

Sims_Result.R: Code to generate simulation result tables

Functions.R: Functions to compute the joint likelihood and Hessian matrix assuming constant baseline hazard function (i.e., exponential model)

Functions_Quadratic.R: Functions to compute the quadratic specification of model needed for Table 4 simulation comparisons. 

Functions_Cubic.R: Functions to compute the cubic specification of model needed for Table 4 simulation comparisons.

Real_Functions_Quadratic.R: Functions to compute the joint profile likelihood and Hessian matrix based on the Breslow approach with left-truncation, for the quadratic epidemic model

Real_Functions_Cubic.R: Functions to compute the joint profile likelihood and Hessian matrix based on the Breslow approach with left-truncation, for the cubic epidemic model

COVID_Quadratic.R: COVID-19 and liver transplant analysis with quadratic epidemic model

COVID_Cubic.R: COVID-19 and liver transplant analysis with cubic epidemic model
