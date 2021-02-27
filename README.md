# qMRI efficiency

This repo contains the code necessary to reproduce the results published in the technical note "Efficiency analysis for quantitative MRI of T1 and T2 relaxometry methods", submitted to PMB (link will be made available to a pre-print).  
The code is divided in several sub-folders, which contain the relevant scripts/functions for each topic. Separately, several auxiliary functions are given in the `/library` sub-folder. Please add this to your Matlab path. Moreover, .mat files are provided as a release; download these and add to your path if you do not want to run all the scripts (some dependencies exist). 
This code is distributed under the MIT licence. If you find it useful please cite the publication once available.

David Leitao, King's College London, January 2021. [@davidmcleitao](https://twitter.com/davidmcleitao)

### Simple experiment

Plots the results from the optimised simple fingerprint experiment with just 5 RF pulses / measurements, as described in the technical note. 

* `graph_simple_MRF_experiments.m` loads the optimised acquisition settings for the simple MRF experiment and plots these. Note that the optimised acquisition settings were obtained using code in the directory `/efficiency_analysis/transient_methods/spoiled_MRF`.

### Efficiency analysis

Contains the code to optimise each analysed method for maximum efficiency and analyses their efficiency distribution for a wide range of tissue parameters - the script `efficiency_comparison.m` performs the comparison of all analysed qMRI methods.  
It is further sub-divided into `/steady-state` and `/transient methods`.

##### Steady-state methods

Contains folders for each steady-state method analysed: `/DESPOT_JSR`, `DESS`, `/PLANET` and `/TESS`.  
Each of these folders contain a script that performs the optimisation of each method acquisition settings for several acquisitions (`optimisation_*.m`), which are saved into a .mat file (`opt_param_*.mat`), and a script that analyses the best acquisition (`analysis_efficiency_*.m`) and saves the results into another .mat file (`*_eff.mat`). Both scripts make use of auxiliary functions provided in the same folder; to calculate the cost function (`cost_function_*.m`) or to help in the analysis of the acquisition settings (`analysis_acq_set_*.m`). 

##### Transient methods

Contains folders for each transient method analysed: `/balanced_MRF` and `/spoiled_MRF`. Each is further sub-divided into `/DE` and `/nonDE`.  
Each of these folders contain two scripts that perform the optimisation of each method acquisition setings for several acquisition: one uses a multi-start strategy (`optimisation_*_multi_start.m`) and the other uses a single-start strategy (`optimisation_*_single_start.m`). The results from each one of them sare saved into .mat files (`opt_param_*_ms.mat` and `opt_param_*_ss.mat`, respectively). The other script analyses the best acquisition (`analysis_effiency_*_.m`) out of the two optimisation results and saves the result into another .mat file (`*_eff.mat`). Analysis and cost function functions were written in .cpp and provided as MEX files in the `/library`.

### Undersampling analysis

Contains the code to perform the approximate dynamics-factor calculation for both spiral and random sampling.  
* `undersampling_impact_random_sampling.m` and `undersampling_impact_spiral_sampling.m` perform the Monte Carlo simulation for random and spiral sampling, respectively, and save the results into several .mat files (`*_sampling_SNR*.mat` and `*_sampling_SaR.mat`, where the latter contains the results of Monte Carlo simulations with no added thermal noise, just aliasing).  
* `spiral_DCF_calculation.m` calculates the DCF for the spiral sampling implemented in `spiral_trajectory.m` and saves the DCF in the .mat file `spiral_DCF.mat`.  
* `analysis_undersampling_impact.m` uses the results from the Monte Carlo simulations to calculate the dynamics-factors.

### Library

The library folder contains several functions used by scripts and functions described above, as well as the [Colormaps Matlab package](https://uk.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps) that provides different colormaps for plotting in Matlab.  
The .cpp code provided in order to be compiled into MEX files requires the [Eigen library v3](http://eigen.tuxfamily.org/index.php?title=Main_Page) to be compiled. Compiled versions for Windows and Linux are already provided.
