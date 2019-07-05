### Instructions for the particle metropolis within Gibbs (PMwG) code in Matlab  ###

* LBA_PMwG_v1.m : This is the primary file to run the algorithm. Start
here. This file loads subsidiary scripts as required, including:
   -> LBA_realdata.mat : Data from Forstmann et al. (2008) in a format ready for estimation with PMwG.
   -> LBA_MC_v1.m : Monte Carlo algorithm for generating initial estimates of the random effects.
   -> LBA_CMC_v1.m : Conditional Monte Carlo algorithm for updating estimates of the random effects.
* LBA_Forstmann_v1.mat : Posterior samples for Forstmann et al. (2008) estimated by PMwG_v1.m.

When applying PMwG to different data sets, the following scripts may
require editing depending on the design of the new data set and the
model that is estimated.
* LBA_n1PDF_reparam_real.m : LBA race equation for two-choice task. 
* LBA_tpdf.m : density function for a single accumulator.
* LBA_tcdf.m : distribution function for a single accumulator.
* reshape_v.m : assign appropriate drift rates (correct, error) to appropriate trials.
* reshape_b.m :  assign appropriate response thresholds (accuracy, neutral, speed) to appropriate trials. 

Subsidiary Matlab functions. These don't require editing for application to different data sets.
* logmvnpdf.m : log density of multivariate normal distribution.
* topdm.m : transform to positive definite symmetric matrix.
* multitransp.m : transposing arrays of matrices.
* multiprod.m : multiplying sub-arrays
