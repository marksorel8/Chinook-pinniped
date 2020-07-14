# Chinook-pinniped

This repository contains code and data used in an analysis of population-specific mortality of endangered Chinook salmon in the Lower Colmbia River associated with increased sea lion abundance in 2013-2015.

The main script for the analysis is "Chinook_sea_lion_mortality_analysis.R". This script conducts the ananlysis but sources several other scripts contained within the "src" folder. Several packages are required to run this script, including "TMB" or "Template Model builder", which also requires RTools. [Instruction here](https://github.com/kaskr/adcomp/wiki/Download). 

The Template MOdel Builder code is contained within "pop_surv_final.cpp".

"data_processing.R" contains functions to process data for the analysis.

"model_fitting_and_results_processing_functions.r" contians functions for fitting models and manipulating outputs.

"plotting_functions.R" contains functions for plotting and reporting analysis results. 
