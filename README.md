# Unbiased-Stationary-MVSDEs
This repository contains the code for the paper titled "Unbiased Approximations for Stationary Distributions of McKean-Vlasov SDEs".  
## How To Use
- The codes `CW.m`, `parameter_estimation.m`, and `n3d.m` produce the estimate introduced in the paper and print it with the line `disp(mean(samples))`.
- To use the other codes, the lines that start with `writematrix(...` in `CW.m`, `parameter_estimation.m`, and `n3d.m` must be uncommented. This will output several files that contain information such as runtime and the particles used to calculate the estimates.
- After these files are output, you can run the other corresponding codes to calculate the estimated density and the MSE rate.
## Comments
While the codes above work on personal computers, the figures in the paper were produced with more than 100K independent samples. To achieve this, the codes were run on the supercomputer IBEX at KAUST. This is also the reason that the variables `job_id` and `proc_id` are defined at the beginning of the codes. These variables are hard-set in the codes here, assuming that the codes will be run on a personal computer. 
