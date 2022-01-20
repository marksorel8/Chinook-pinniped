# Chinook-pinniped

This repository contains code and data used in the manuscript:

[Sorel, Mark H., Richard W. Zabel, Devin S. Johnson, A. Michelle Wargo Rub, and Sarah J. Converse. "Estimating population‐specific predation effects on Chinook salmon via data integration." Journal of Applied Ecology 58, no. 2 (2021): 372-381.](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/1365-2664.13772)


The main script for the analysis is `Chinook_sea_lion_mortality_analysis.R`. This script conducts the ananlysis but sources several other scripts contained within the `src` folder. Several packages are required to run this script, including "TMB" (Template Model builder), which  requires RTools. [Download instruction here](https://github.com/kaskr/adcomp/wiki/Download). 

The following files, which contian code for the anlaysis, can be found in the `src` folder.

`pop_surv_final.cpp` contains the Template Model Builder code for the HMM and the population-specfic migration timing.

`data_processing.R` contains functions to process data for the analysis.

`model_fitting_and_results_processing_functions.r` contians functions for fitting models and manipulating outputs.

`plotting_functions.R` contains functions for plotting and reporting analysis results. 

If you have questions, please contact me through this repository or at marks6@uw.edu. Thank you. 
