# PRO-ADD
This repository contains the code used to simulate results for the paper **PRO-ADD: Patient-Empowered Dose-Finding Trials by Integrating Safety, Efficacy and Patient-Reported Outcomes for Optimal Dose Selection**. 

## Background 
This project evaluates the operating characteristics of the novel PRO-ADD design under simulation scenarios.

## Description of R files
* **functions.R** - code for Section 4: functions required to generate trials and identify the OBD for 5,000 simulations running in parallel. 
  
* **sim.R** - code for Section 4: code required to run PRO-ADD simulation studies. Functions required to run this code are defined in `functions.R`.

* **shape_param_inc.csv** - code for Section 4: .csv file containing matrix of shape parameters to define PRO-nAE burden score simualation scenarios.

* **rate_param_inc.csv** - code for Section 4: .csv file containing rate parameter to define PRO-nAE burden score simualation scenarios.
