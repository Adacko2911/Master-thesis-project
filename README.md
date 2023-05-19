## Should you make the Explicit-duration Hidden Markov model your asset when analysing intensive longitudinal data for behavioral research? 
#### Methodology and Statistics for the Behavioural, Biomedical and Social Sciences (MSBBSS), Utrecht University
### Aleksandra Dacko 

#### Reasearch archive
 
---
This repository contains R code to reproduce the results of the master thesis. 

### Repository structure
In the following research archive you can find the following folders[^1]:
```
MHMM_vs_MEDHMM/
├── README.md
├── Session_info.txt - includes information regarding all packages loaded duering the analysis as well as the R software specifications
├── data/
│   ├── data_simulation.R - simulated data generating script
│   └── read_me.txt  - description of the data files that should be included
├── convergence_run/   
│    ├── medhmm/
│    │   ├── convergence_run_medhmm*.R - main scripts running the convergence inspections
│    │   ├── doit* - doit files necessary to run to submit jobs included in jobs/ folder
│    │   ├── jobs/ - jobs that have been submitted 
│    │   ├── outputs/  - all outputs are stored in the folder including pdfs of the runs that had been included ad hoc for reader inspection
│    │   └── start_convergence_medhmm*.sh - .sh files used to submit jobs to snellius
│    │   
│    └── mhmm/
│        ├── convergence_runs.R - main scripts running the convergence inspections
│        ├── doit - doit files necessary to run to submit jobs included in jobs/ folder
│        ├── jobs/ - jobs that have been submitted 
│        ├── outputs/  - all outputs are stored in the folder including pdfs of the runs that had been included ad hoc for reader inspection
│        └── start_convergence_mhmm.sh - .sh file used to submit jobs to snellius
├── medhmm_run/
│    ├── model_fitting2_and_post_processes.R - main scripts running the simulation study for the multilevel explicit duration hidden markov model
│    ├── model_fitting2_and_post_processes_extra_run.R - main scripts running the additional scenario of simulation study for the multilevel explicit duration hidden markov model
│    ├── doit* - doit files necessary to run to submit jobs included in jobs/ folder
│    ├── jobs/ - jobs that have been submitted
│    ├── outputs/  - all outputs are stored in the folder including including data that is later used to obtain the statistics
│    ├── start_2dep_medhmm_*.sh - .sh files used to submit jobs to snellius
│    └── parameters/ - folder including all parameters of the simulation scenarios that are further passed to snellius engine
│
├── mhmm_run/
│    ├── model_fitting2_and_post_processes.R - main scripts running the simulation study for the multilevel  hidden markov model
│    ├── model_fitting2_and_post_processes_extra.R - main scripts running the additional scenario of simulation study for the multilevel hidden markov model
│    ├── doit* - doit files necessary to run to submit jobs included in jobs/ folder
│    ├── jobs/ - jobs that have been submitted
│    ├── outputs/  - all outputs are stored in the folder including including data that is later used to obtain the statistics
│    ├── start_2dep_medhmm_*.sh - .sh files used to submit jobs to snellius
│    └── parameters/ - folder including all parameters of the simulation scenarios that are further passed to snellius engine
│
├── simulation_main/
│    ├── parse outputs.R - script parsing outputs included in mhmm_run/outputs/ and medhmm_run/outputs/ that saves the final post processed files in simulation_main/post_processed_data/
│    ├── plot_scripts/ - includes three files that generate the plots that are further used in the manuscript 
│    ├── outputs plots/ - includes all plots generate by the scripts included in simulation_main/plot_scripts/
│    └── post_processed_data/ include final outputs that are further used for table generating and plot building. Main Results.xlsx inluded in this folder has been manually parsed and included for and easy overview.
│
├── manuscript/ - includes all laTeX files and graphics that were used to generate the manuscript of the thesis
│
└── example/ 
    ├── outputs/ - includes the outputs of the empirical example model fitting 
    ├── plots/ - includes all plots included in the final manuscript that are 
    ├── train_models.R - code to train both MEDHMM and MHMM with a use of the empirical data
    ├── models_postprocess.R - scripts generating trace plots and other statistics reported in main manuscript. In addition includes code to generate figures 9, 10, 12, C4. Also includes script to simulate data for the posterior predictive checks 
    ├── summarise_data.R - scripts summarizing in sample data to obtain the group-level and patient-level statistics. The data is further used for the posterior predictive checks 
    ├── utility_functions.R - script including function to summaries, simulate and plot data
    ├── ppc_mean_code.R - code to generate the plots of the posterior predictive checks i.e. figures C1 and C2 from the manuscript
    └── decoding_figures_11_C5.R - script generation figures 11 and C5 from the manuscript

 ```   
[^1]: The "*" symbol indicates a wildcard for the file names as a lot of them were generated with the same structure. 

## Reproducing Results

To reproduce the results of the study, follow these steps:

1. Download the code in the `convergence_run`,`example`,`medhmm_run`,`mhmm_run` and `data` folder and the empirical data (upon request from Groningen Medical centre).
2. Run the R scripts in the following order to reproduce the results from main simulation study:
    1. `data/data_simulation.R`: simulate the data.
    2. `convergence_run/medhmm/convergence_run_medhmm1.R`,`convergence_run/medhmm/convergence_run_medhmm2.R`, `convergence_run/mhmm/convergence_runs.R` : analyse the convergence results from fitting the Bayesian multilevel EDHMM and HMM on the simulated data.
    3. R (.R) scripts included in `medhmm_run/` and `mhmm_run/` directories: run the simulation.
    4. `simulatipon_main/parse_outputs.R`: summarize and parse all results. (results `produced to ./post_processed_data/`).
    5. R (.R) scripts included in `simulatipon_main/plot_scripts/`directory: produce figures (results `./produced to plot_scripts/`).
3. Run the R scripts in the following order to reproduce the results from empirical exmple:
    1. `example/train_models.R`: fit the data to MEDHMM and the MHMM.
    2. `example/summarise_data.R`: get the summary statistics of individuals and study group.
    3. `example/utility_functions.R`: load functions needed.
    4. `example/models_postprocess.R`: summarize data results, produce plots and convergence trace plots.
    5. `example/decoding_figures_11_C5.R` and `example/ppc_mean_code.R`: to obtain figures and posterior predictive check analisys results.
  
Notice that the simulation is computationally intensive; it was run in the [Dutch National Supercomputer Snellius](https://www.surf.nl/en/dutch-national-supercomputer-snellius) using a single node with varying number of cores. The cluster computer scripts are available in `convergence_run/medhmm/`, `convergence_run/mhmm/`,`medhmm_run/`, and `mhmm_run/` and consists of the **job** files, **doit** files and executable **.sh** files. Those files can be reuse however they need to be adjusted in order to mach future user file directories. 



## Contact

The empirical data cannot be found in this repository due to privacy issues. However, the data can be made available upon reasonable request from researchers. Please contact <f.m.bos01@umcg.nl> for more information.

For any further questions, please contact <a.dacko@students.uu.nl>.
