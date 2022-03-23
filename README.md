# Fed-GLMM

Fed-GLMM enables the joint implementation of generalized Linear mixed models (GLMMs) for datasets from multiple sites. It requires three steps:

1. Fit GLMM locally to obtain initial values of parameter estimates
2. Calculate and broadcast the summary statistics which are the local first- and second-order derivatives
3. Update parameter estimates using the summary statistics

When iterative communications are allowed, steps 2-3 can then be repeated to further update the parameter values. 
