# Manuscript code repository
This repository contains statistical code related to peer-reviewed published manuscripts. The goal is to increase transparency and reproducibility in research.

## Zabor_Validity-of-a-method-for-etiologic-heterogeneity

Files:

- `eh_cluster_sim_fun.R` contains the primary function for conducting a clustering simulation
- `run_eh_cluster_sim_fun.R` is where you can set the simulation parameters to desired values and run the simulation

Instructions:

1. Download both files to the same directory
2. Change the filepath in `setwd()` to point to the directory where the files were downloaded
3. Set simulation parameters in `run_eh_cluster_sim.R` to desired values and run
4. Use `set.seed()` to set a seed. Note that a variety of seeds were used to produce results in the manuscript.
