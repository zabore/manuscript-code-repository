# Manuscript code repository

This repository contains statistical code and data related to published manuscripts and book chapters. The goal is to increase transparency and reproducibility in research.

***

## Zabor-Eaton_ipcw-tutorial

Reference:
>

Code for the guided example and the simulation studies, for both single and competing events.

File:

- `single-event-guided-example.R` code for the single event guided example data and analysis
- `single_example_dat.rda` single simulated dataset for the single event guided example
- `single_example_ipcw_dat.rda` single simulated dataset for the single event guided example, in long format with ipc weights
- `single-event-functions.R` functions used in both the single event guided example and simulation study
- `single-event-simulation-study.R` code to generate simulation study results for single event setting
- - `simulation_code_cuminc_finegray.R` code to generate simulation study results for the competing events setting

***

## Mi-Zabor_bootstrap-impute-predict-tutorial

Reference:
>

Code to simulate and analyze the single synthetic dataset used in the guided example.

File:

- `sim-guided-example-data.R` code to generate single synthetic dataset with missing values used in guided example
- `dat0.rda` single synthetic dataset with missing values used in guided example, generated by `sim-guided-exmaple-data.R` and used in `analyze-guided-example.R`
- `analyze-guided-example.R` code to bootstrap, impute, and obtain performance measures on the single synthetic dataset used in the guided example
- `simulation-functions.R` collection of functions to run a simulation study
- `simulation-studies.R` code to run simulation studies included in the paper

***

## Zabor_Clinical-Opthalmic-Oncology-Cancer-Survival

Reference:
> Clinical Ophthalmic Oncology: Basic Principles. United States, Springer.

Synthetic data file, R code, and resulting Quarto report associated with all results included in the book chapter are available:

Files:

- `uveal_survival_data.csv` synthetic dataset in csv format
- `uveal_survival_data.rds` synthetic dataset in rds format
- `cancer_survival_report.qmd` R Quarto file containing R code to produce all results presented in the chapter
- `cancer_survival_report.docx` Word report rendered from the associated .qmd file

***

## Zabor_Randomized-Biomarker-Guided Designs

Reference:
> Zabor EC, Kaizer AM, Pennell NA, Hobbs BP. Optimal predictive probability designs for randomized biomarker-guided oncology trials. Front Oncol. 2022;12:955056. Published 2022 Dec 6. doi:10.3389/fonc.2022.955056

Data, data generating scripts, and analysis scripts for the three designs are available:

1. Pooled control arm design
    - `data\1-pooled-generate-data.R` script to generate the simulated trial data for the pooled control arm design
    - `data\p-sim-dat.rda` resulting simulated pooled control arm trial data
    - `analysis\2-pooled-apply-ppseq.R` script to run the analysis for the pooled control arm design
    - `analysis\p-ppseq-res.rda` analysis results to obtain operating characteristics for the pooled control arm design
2. Stratified control arm design
    - `data\1-stratified-generate-data.R` script to generate the simulated trial data for the stratified control arm design
    - `data\s-sim-dat.rda` resulting simulated stratified control arm trial data
    - `analysis\2-stratified-apply-ppseq.R` script to run the analysis for the stratified control arm design
    - `analysis\s-ppseq-res.rda` analysis results to obtain operating characteristics for the stratified control arm design
3. Enrichment design
    - `analysis\1-enrichment-apply-ppseq.R` script to run the analysis for the enrichment design
    - `analysis\e-ppseq-res.rda` analysis results to obtain operating characteristics for the enrichment design

Notes:
- these files were organized in an elaborate folder structure with a top-level R project and all filepaths used `here::here` relative to the R project. Please alter the filepaths accordingly
- no data file or data generating script is included with the enrichment design since data from the pooled control group design were utilized for stage 1 and data from the stratified control group design were utilized for stage 2

***

## Eaton-Zabor_Interval-RFS-sim-study

Reference:
> Eaton AA, Zabor EC. Analysis of composite endpoints with component-wise censoring in the presence of differential visit schedules. Stat Med. 2022 Apr 30;41(9):1599-1612. doi: 10.1002/sim.9312. Epub 2022 Jan 18. PMID: 35043427.

Files:

- `fn-create-single-dataset.R` function to simulate a single interval censored dataset
- `fn-fit-models.R` function to fit each of the models included in the paper
- `fn-tidy-ic_sp.R` helper function to tidy the results from icenReg::ic_sp() function
- `program1-run-simulation.R` program to run a simulation for a given scenario and summarize results

Instructions:
1. Download all files to the same directory
2. Open `program1-run-simulation.R`
3. Install the `cwcens` package in R from `remotes::install_github("anneae/cwcens")`
4. Either set your working directory or alter the file paths for the 3 function files in the section "Load functions"
5. Change the parameters of interest in the section "Generate many datasets", including possibly the number of simulated datasets and the seed
6. In the section "Prepare data to summarize" you will need to change the truth to appropriately match the setting you are examining. In the example the HR for both recurrence and death were set to 1.5 so the true log HR in all cases is log(1.5)
7. Run all code in the file to fit the models and summarize the results with a boxplot and table

***

## Zabor_Basket-Trial-with-FDR-Control

Reference:
> Zabor EC, Kane MJ, Roychoudhury S, Nie L, Hobbs BP. Bayesian basket trial design with false-discovery rate control. Clin Trials. 2022 Feb 7:17407745211073624. doi: 10.1177/17407745211073624. Epub ahead of print. PMID: 35128970.

This paper contains three distinct parts:

1. Model calibration
    - `heterogeneous_response_sim_fn.R` A function to run the basket trial with different design parameters
    - `heterogeneous_response_sim_do.R` A function that implements `heterogeneous_response_sim_fn.R` for the design parameters of interest
    - `heterogeneous_response_sim_res.rda` The simulation results from `heterogeneous_response_sim_do.R`
    - `heterogeneous_response_sim_report.Rmd` R Markdown report to summarize the results in `heterogeneous_response_sim_res.rda`
    
2. Case study
    - `neratinib_mem_test.R` The program to generate simulation results based on the neratinib case study
    - `nerat_basket_0.25.rda` and `nerat_basket_0.5.rda` are two files of simulation results generated by `neratinib_mem_test.R` for different prior probabilities of exchangeability
    - `neratinib_mem_results.Rmd` R Markdown report to summarize the results in `nerat_basket_0.25.rda` and `nerat_basket_0.5.rda`
    
3. Trial operating characteristics
    - `do_single_sim.R` A function to generate MEM and frequentist results for baskets with given sample sizes and true response rates
    - `scenario1-alt0.R` Executes the global alternative scenario
    - `scenario1-alt0.rda` Results from `scenario1-alt0.R`
    - `scenario1-alt1.R` Executes alternative 1 scenario
    - `scenario1-alt1.rda` Results from `scenario1-alt1.R`
    - `scenario1-alt2.R` Executes alternative 2 scenario
    - `scenario1-alt2.rda` Results from `scenario1-alt2.R`
    - `scenario1-alt3.R` Executes alternative 3 scenario
    - `scenario1-alt3.rda` Results from `scenario1-alt3.R`
    - `scenario1-null.R` Executes the global null scenario
    - `scenario1-null.rda` Results from `scenario1-null.R`
    - `operating-chars-results.R` Gets summary dataframes of the results from the `scenario1-****.rda` files
    - `scenario1-summary.rda` Summary results from `operating-chars-results.R`
    - `operating-chars-results-report.Rmd` R Markdown report to summarize the results in `scenario1-summary.rda`
    
*Please note that these files were organized in an elaborate folder structure with a top-level R project and all filepaths used `here::here` relative to the R project. Please alter the filepaths accordingly.*

***

## Zabor_Validity-of-a-method-for-etiologic-heterogeneity

Reference:
> Zabor EC, Seshan VE, Wang S, Begg CB. Validity of a method for identifying disease subtypes that are etiologically heterogeneous. Stat Methods Med Res. 2021 Sep;30(9):2045-2056. doi: 10.1177/09622802211032704. Epub 2021 Jul 28. PMID: 34319833.

Files:

- `eh_cluster_sim_fun.R` contains the primary function for conducting a clustering simulation
- `run_eh_cluster_sim_fun.R` is where you can set the simulation parameters to desired values and run the simulation

Instructions:

1. Download both files to the same directory
2. Install the `riskclustr` package from CRAN using `install.packages("riskclustr")`. Note that this line of code is commented out in `eh_cluster_sim_fun.R` as it only needs to be run once.
3. Change the filepath in `setwd()` in `run_eh_cluster_sim.R`to point to the directory where the files were downloaded
4. Set simulation parameters in `run_eh_cluster_sim.R` to desired values 
5. Use `set.seed()` to set a seed. Note that a variety of seeds were used to produce results in the manuscript.
6. Run `run_eh_cluster_sim.R` to produce results. See notes at bottom of code file for details on results structure.
